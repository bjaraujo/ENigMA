// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <iostream>
#include <algorithm>

#include "GeoCircle.hpp"
#include "GeoLine.hpp"
#include "MshMeshQuery.hpp"
#include "MshTriangle.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        CMshTriangleMesher<Real>::CMshTriangleMesher()
            : m_begin()
            , m_end()
            , m_nextEdgeId(0)
            , m_previousNbElements(0)
            , m_timeInterval(1)
            , m_dataInterval(100)
            , m_bStop(false)
            , onUpdate(nullptr)
        {
        }

        template <typename Real>
        CMshTriangleMesher<Real>::~CMshTriangleMesher()
        {
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::checkUpdate()
        {
            if (onUpdate != nullptr)
            {
                m_end = clock();

                if (static_cast<Real>(m_end - m_begin) / CLOCKS_PER_SEC > m_timeInterval)
                {
                    m_begin = m_end;

                    if (m_surfaceMesh.nbElements() > m_previousNbElements + m_dataInterval)
                    {
                        m_previousNbElements = m_surfaceMesh.nbElements();
                        onUpdate(true);
                    }
                    else
                        onUpdate(false);
                }
            }
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::removeEdge(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance)
        {
            anAdvEdge.remove = true;
            removeEdgeFromRtree(anAdvEdge, aTolerance);
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::addEdgeToRtree(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance)
        {
            Integer aNodeId1 = anAdvEdge.nodeId[0];
            Integer aNodeId2 = anAdvEdge.nodeId[1];

            CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
            CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);

            CGeoBoundingBox<Real> aBoundingBox;

            aBoundingBox.addCoordinate(aNode1);
            aBoundingBox.addCoordinate(aNode2);

            aBoundingBox.grow(aTolerance);

            m_tree.addGeometricObject(anAdvEdge.id, aBoundingBox);
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::removeEdgeFromRtree(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance)
        {
            Integer aNodeId1 = anAdvEdge.nodeId[0];
            Integer aNodeId2 = anAdvEdge.nodeId[1];

            CGeoBoundingBox<Real> aBoundingBox;

            aBoundingBox.addCoordinate(m_surfaceMesh.node(aNodeId1));
            aBoundingBox.addCoordinate(m_surfaceMesh.node(aNodeId2));

            aBoundingBox.grow(aTolerance);

            m_tree.removeGeometricObject(anAdvEdge.id, aBoundingBox);
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::edgeExists(SMshAdvancingFrontEdge<Real>& anAdvEdge, Integer& aDuplicateEdgeId, std::vector<Integer>& sEdges)
        {
            std::vector<Integer> sNodeIds;

            sNodeIds.push_back(anAdvEdge.nodeId[0]);
            sNodeIds.push_back(anAdvEdge.nodeId[1]);

            std::sort(sNodeIds.begin(), sNodeIds.end());

            for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j)
            {
                // Discard same edge
                if (sEdges[j] == anAdvEdge.id)
                    continue;

                SMshAdvancingFrontEdge<Real>& anotherEdge = m_anAdvFront[sEdges[j]];

                std::vector<Integer> sOtherNodeIds;

                sOtherNodeIds.push_back(anotherEdge.nodeId[0]);
                sOtherNodeIds.push_back(anotherEdge.nodeId[1]);

                std::sort(sOtherNodeIds.begin(), sOtherNodeIds.end());

                if (sOtherNodeIds[0] == sNodeIds[0] && sOtherNodeIds[1] == sNodeIds[1])
                {
                    aDuplicateEdgeId = sEdges[j];
                    return true;
                }
            }

            return false;
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::edgeOk(SMshAdvancingFrontEdge<Real>& anAdvEdge, CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, std::vector<Integer>& sEdges, const Real aTolerance)
        {
            CGeoCoordinate<Real> aPoint;

            CGeoLine<Real> aLine1(aNode1, aNode2);

            for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j)
            {
                // Exclude current edge
                if (sEdges[j] == anAdvEdge.id)
                    continue;

                SMshAdvancingFrontEdge<Real>& anotherEdge = m_anAdvFront[sEdges[j]];

                if (anotherEdge.remove)
                    continue;

                CGeoLine<Real>& aLine2 = anotherEdge.line;

                CGeoIntersectionType anIntersectionType;

                if (aLine1.intersects(aLine2, aPoint, anIntersectionType, aTolerance))
                {
                    if (anIntersectionType == IT_INTERNAL)
                        return false;
                }
            }

            return true;
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::findClosestNodes(std::vector<Integer>& sEdges, std::vector<Integer>& sNodes)
        {
            sNodes.clear();

            for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j)
            {
                if (m_anAdvFront[sEdges[j]].remove)
                    continue;

                for (Integer k = 0; k < 2; ++k)
                {
                    Integer aNodeId = m_anAdvFront[sEdges[j]].nodeId[k];

                    if (std::find(sNodes.begin(), sNodes.end(), aNodeId) == sNodes.end())
                        sNodes.push_back(aNodeId);
                }
            }

            for (Integer j = 0; j < static_cast<Integer>(m_interiorNodes.size()); ++j)
            {
                if (m_interiorNodes[j].remove)
                    continue;

                Integer aNodeId = m_interiorNodes[j].nodeId;

                if (std::find(sNodes.begin(), sNodes.end(), aNodeId) == sNodes.end())
                    sNodes.push_back(aNodeId);
            }
        }

        template <typename Real>
        Real CMshTriangleMesher<Real>::findShortestDistance(std::vector<Integer>& sEdges, CGeoLine<Real>& aLine, Integer anAdvEdgeId, const Real aTolerance)
        {
            Real dmin = std::numeric_limits<Real>::max();

            for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j)
            {
                if (sEdges[j] == anAdvEdgeId)
                    continue;

                if (m_anAdvFront[sEdges[j]].remove)
                    continue;

                CGeoCoordinate<Real> aNewPoint;

                CGeoIntersectionType anIntersectionType;

                if (m_anAdvFront[sEdges[j]].line.intersects(aLine, aNewPoint, anIntersectionType, aTolerance))
                {
                    Real distance = (aNewPoint - aLine.startPoint()).norm();

                    if (distance < dmin)
                        dmin = distance;
                }
            }

            return dmin;
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::triangleContainsNode(CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance)
        {
            CMshTriangle<Real> aTriangle;

            aTriangle.addVertex(aNode1);
            aTriangle.addVertex(aNode2);
            aTriangle.addVertex(aNode3);

            for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j)
            {
                aNodeId = sNodes[j];

                CMshNode<Real>& aNode = m_surfaceMesh.node(aNodeId);

                if ((aNode - aNode1).norm() < aTolerance || (aNode - aNode2).norm() < aTolerance || (aNode - aNode3).norm() < aTolerance)
                    continue;

                CGeoIntersectionType anIntersectionType;

                if (aTriangle.contains(aNode, anIntersectionType, aTolerance))
                {
                    if (anIntersectionType == IT_INTERNAL || anIntersectionType == IT_EDGE)
                        return true;
                }
            }

            return false;
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::checkDelaunay(CMshNode<Real>& aNewNode, const Real aTolerance)
        {
            for (Integer i = 0; i < m_surfaceMesh.nbElements(); ++i)
            {
                Integer anElementId = m_surfaceMesh.elementId(i);

                CMshElement<Real>& anElement = m_surfaceMesh.element(anElementId);

                if (anElement.elementType() != ET_TRIANGLE)
                    continue;

                Integer aNodeId1 = anElement.nodeId(0);
                Integer aNodeId2 = anElement.nodeId(1);
                Integer aNodeId3 = anElement.nodeId(2);

                CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
                CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);
                CMshNode<Real>& aNode3 = m_surfaceMesh.node(aNodeId3);

                CGeoCircle<Real> aCircle(aNode1, aNode2, aNode3);

                if (aCircle.contains(aNewNode, aTolerance))
                    return false;
            }

            return true;
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::adjustConnectivity(std::vector<Integer>& sEdges)
        {
            for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i)
            {
                Integer anAdvEdgeId = i;

                SMshAdvancingFrontEdge<Real>& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                if (anAdvEdge.remove)
                    continue;

                if (anAdvEdge.nodeId[0] == anAdvEdge.nodeId[1])
                {
                    anAdvEdge.remove = true;
                    continue;
                }

                if (anAdvEdge.neighborId[0] >= static_cast<Integer>(m_anAdvFront.size()) || anAdvEdge.neighborId[1] >= static_cast<Integer>(m_anAdvFront.size()))
                {
                    throw std::out_of_range("Connectivity is out of range!");
                }

                SMshAdvancingFrontEdge<Real>& aPrevEdge = m_anAdvFront[anAdvEdge.neighborId[0]];
                SMshAdvancingFrontEdge<Real>& aNextEdge = m_anAdvFront[anAdvEdge.neighborId[1]];

                if (aPrevEdge.remove || aNextEdge.remove)
                {
                    for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j)
                    {
                        Integer anotherEdgeId = sEdges[j];

                        // Discard same edge
                        if (anotherEdgeId == anAdvEdge.id)
                            continue;

                        SMshAdvancingFrontEdge<Real>& anotherEdge = m_anAdvFront[anotherEdgeId];

                        if (anotherEdge.remove)
                            continue;

                        if (anAdvEdge.nodeId[0] == anotherEdge.nodeId[1])
                        {
                            anAdvEdge.neighborId[0] = anotherEdgeId;
                            anotherEdge.neighborId[1] = anAdvEdgeId;
                        }
                        else if (anAdvEdge.nodeId[1] == anotherEdge.nodeId[0])
                        {
                            anAdvEdge.neighborId[1] = anotherEdgeId;
                            anotherEdge.neighborId[0] = anAdvEdgeId;
                        }
                    }
                }
            }
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::cleanDuplicateEdges(std::vector<Integer>& sEdges, const Real aTolerance)
        {
            if (m_nextEdgeId < 2)
                return;

            bool bAdjustConnectivity = false;

            // Check last two edges
            for (Integer i = 0; i < 2; ++i)
            {
                Integer anAdvEdgeId = m_nextEdgeId - i - 1;

                SMshAdvancingFrontEdge<Real>& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                if (anAdvEdge.remove)
                    continue;

                Integer aDuplicateEdgeId;

                if (edgeExists(anAdvEdge, aDuplicateEdgeId, sEdges))
                {
                    SMshAdvancingFrontEdge<Real>& aDuplicateEdge = m_anAdvFront[aDuplicateEdgeId];

                    this->removeEdge(anAdvEdge, aTolerance);
                    this->removeEdge(aDuplicateEdge, aTolerance);

                    bAdjustConnectivity = true;
                }
            }

            if (bAdjustConnectivity)
                this->adjustConnectivity(sEdges);
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::addTriangle(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Integer aNodeId, std::vector<Integer>& sEdges, const Real aTolerance)
        {
            Integer aNodeId1 = anAdvEdge.nodeId[0];
            Integer aNodeId2 = anAdvEdge.nodeId[1];
            Integer aNodeId3 = aNodeId;

            if (aNodeId3 == aNodeId1 || aNodeId3 == aNodeId2)
            {
                std::cout << "Error: invalid triangle!" << std::endl;
                return;
            }

            // Add new triangle
            CMshElement<Real> aNewElement(ET_TRIANGLE);
            aNewElement.addNodeId(aNodeId1);
            aNewElement.addNodeId(aNodeId2);
            aNewElement.addNodeId(aNodeId3);

            Integer aNewElementId = m_surfaceMesh.nextElementId();

            m_surfaceMesh.addElement(aNewElementId, aNewElement);

            SMshAdvancingFrontEdge<Real>& aPrevEdge = m_anAdvFront[anAdvEdge.neighborId[0]];
            SMshAdvancingFrontEdge<Real>& aNextEdge = m_anAdvFront[anAdvEdge.neighborId[1]];

            Integer aNewEdgeId1 = m_nextEdgeId++;
            Integer aNewEdgeId2 = m_nextEdgeId++;

            // Add edge 1
            SMshAdvancingFrontEdge<Real> aNewEdge1;
            aNewEdge1.id = aNewEdgeId1;
            aNewEdge1.remove = false;
            aNewEdge1.boundary = false;
            aNewEdge1.nodeId[0] = aNodeId1;
            aNewEdge1.nodeId[1] = aNodeId3;
            aNewEdge1.neighborId[0] = anAdvEdge.neighborId[0];
            aNewEdge1.neighborId[1] = aNewEdgeId2;
            aNewEdge1.elementId = aNewElementId;
            aNewEdge1.nodeNotId3 = aNodeId2;

            aNewEdge1.build(m_surfaceMesh);

            // Correct connectivity
            aPrevEdge.neighborId[1] = aNewEdgeId1;

            // Add edge to rtree
            addEdgeToRtree(aNewEdge1, aTolerance);

            // Add edge 2
            SMshAdvancingFrontEdge<Real> aNewEdge2;
            aNewEdge2.id = aNewEdgeId2;
            aNewEdge2.remove = false;
            aNewEdge2.boundary = false;
            aNewEdge2.nodeId[0] = aNodeId3;
            aNewEdge2.nodeId[1] = aNodeId2;
            aNewEdge2.neighborId[0] = aNewEdgeId1;
            aNewEdge2.neighborId[1] = anAdvEdge.neighborId[1];
            aNewEdge2.elementId = aNewElementId;
            aNewEdge2.nodeNotId3 = aNodeId1;

            aNewEdge2.build(m_surfaceMesh);

            // Correct connectivity
            aNextEdge.neighborId[0] = aNewEdgeId2;

            // Add edge to rtree
            this->addEdgeToRtree(aNewEdge2, aTolerance);

            // Remove this edge
            this->removeEdge(anAdvEdge, aTolerance);

            m_anAdvFront.push_back(aNewEdge1);
            m_anAdvFront.push_back(aNewEdge2);

            this->cleanDuplicateEdges(sEdges, aTolerance);
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::remesh(ENigMA::mesh::CMshMesh<Real>& aMesh, Real meshSize, const Real aTolerance)
        {
            std::stringstream ss(std::stringstream::in | std::stringstream::out);

            ss << meshSize;

            ENigMA::analytical::CAnaFunction<Real> aAnaFunction;

            aAnaFunction.set(ss.str());

            return this->remesh(aMesh, aAnaFunction, aTolerance);
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::remesh(ENigMA::mesh::CMshMesh<Real>& aMesh, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, const Real aTolerance)
        {
            Real x, y;

            meshSizeFunc.removeAllVariables();

            meshSizeFunc.defineVariable("x", x);
            meshSizeFunc.defineVariable("y", y);

            std::vector<bool> processedElements;
            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {
                processedElements.push_back(false);
            }

            // Split edge mesh according to local mesh size
            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {
                Integer anElementId = aMesh.elementId(i);
                CMshElement<Real>& anElement = aMesh.element(anElementId);
                
                if (anElement.elementType() == EElementType::ET_BEAM)
                {
                    Integer aNodeId1 = anElement.nodeId(0);
                    Integer aNodeId2 = anElement.nodeId(1);

                    CMshNode<Real> aNode1 = aMesh.node(aNodeId1);
                    CMshNode<Real> aNode2 = aMesh.node(aNodeId2);

                    bool bConnectivity = true;

                    Integer aFaceId1 = 0;
                    Integer aFaceId2 = 0;

                    if (anElement.nbFaceIds() == 2)
                    {
                        aFaceId1 = anElement.faceId(0);
                        aFaceId2 = anElement.faceId(1);
                    }
                    else
                        bConnectivity = false;

                    CGeoVector<Real> v = aNode2 - aNode1;

                    CMshNode<Real> aMidNode = (aNode1 + aNode2) * 0.5;

                    x = aMidNode.x();
                    y = aMidNode.y();

                    Real localMeshSize = meshSizeFunc.evaluate();

                    Integer ne = static_cast<Integer>(floor(v.norm() / localMeshSize + 0.5));

                    if (ne > 0)
                    {
                        Real de = v.norm() / ne;

                        CGeoNormal<Real> n = v;

                        n.normalize();

                        Integer aPrevNodeId = aNodeId1;
                        Integer aPrevFaceId = aFaceId1;

                        for (Integer j = 0; j < ne - 1; ++j)
                        {
                            // Create new node
                            CMshNode<Real> aNode = aNode1 + (j + 1) * de * n;

                            Integer aNewNodeId = aMesh.nextNodeId();
                            aMesh.addNode(aNewNodeId, aNode);

                            CMshElement<Real> aNewElement(ET_BEAM);
                            aNewElement.addNodeId(aPrevNodeId);
                            aNewElement.addNodeId(aNewNodeId);

                            Integer aNewElementId = aMesh.nextElementId();

                            if (bConnectivity)
                            {
                                std::vector<CMshFace<Real>> sFaces;
                                aNewElement.generateFaces(sFaces);

                                for (Integer k = 0; k < static_cast<Integer>(sFaces.size()); ++k)
                                {
                                    CMshFace<Real> aNewFace = sFaces[k];

                                    Integer aNewFaceId = aMesh.nextFaceId();

                                    aNewFace.setElementId(aNewElementId);

                                    if (k == 0)
                                        aNewFace.setPairFaceId(aPrevFaceId);
                                    else
                                        aNewFace.setPairFaceId(aNewFaceId + 1);

                                    aMesh.addFace(aNewFaceId, aNewFace);
                                    aNewElement.addFaceId(aNewFaceId);

                                    aPrevFaceId = aNewFaceId;
                                }
                            }

                            aPrevNodeId = aNewNodeId;

                            aMesh.addElement(aNewElementId, aNewElement);

                            // Split element
                            aMesh.element(anElementId).setNodeId(0, aNewNodeId);
                            aMesh.element(anElementId).setNodeId(1, aNodeId2);

                            if (bConnectivity)
                            {
                                aMesh.element(anElementId).setFaceId(0, aPrevFaceId);
                                aMesh.element(anElementId).setFaceId(1, aFaceId2);
                            }
                        }
                    }
                }
                else if (anElement.elementType() == EElementType::ET_TRIANGLE)
                {
                    for (Integer j = 0; j < anElement.nbFaceIds(); ++j)
                    {
                        Integer aFaceId = anElement.faceId(j);
                        CMshFace<Real>& aFace = aMesh.face(aFaceId);

                        if (aFace.hasPair())
                        {
                            CMshFace<Real>& aPairFace = aMesh.face(aFace.pairFaceId());
                            Integer aPairElementId = aPairFace.elementId();
                            CMshElement<Real>& aPairElement = aMesh.element(aPairElementId);

                            Integer aNodeId1 = aFace.nodeId(0);
                            Integer aNodeId2 = aFace.nodeId(1);

                            CMshNode<Real> aNode1 = aMesh.node(aNodeId1);
                            CMshNode<Real> aNode2 = aMesh.node(aNodeId2);

                            CGeoVector<Real> v = aNode2 - aNode1;
                            CMshNode<Real> aMidNode = (aNode1 + aNode2) * 0.5;

                            x = aMidNode.x();
                            y = aMidNode.y();

                            Real localMeshSize = meshSizeFunc.evaluate();

                            if (v.norm() > localMeshSize)
                            {
                                Integer aNodeId3 = -1;
                                for (int k = 0; k < anElement.nbNodeIds(); k++)
                                {
                                    if (anElement.nodeId(k) != aNodeId1 && anElement.nodeId(k) != aNodeId2)
                                    {
                                        aNodeId3 = anElement.nodeId(k);
                                        break;
                                    }
                                }

                                Integer aNodeId4 = -1;
                                for (int k = 0; k < aPairElement.nbNodeIds(); k++)
                                {
                                    if (aPairElement.nodeId(k) != aNodeId1 && aPairElement.nodeId(k) != aNodeId2)
                                    {
                                        aNodeId4 = aPairElement.nodeId(k);
                                        break;
                                    }
                                }

                                if (aNodeId3 != -1 && aNodeId4 != -1 && !processedElements[anElementId] && !processedElements[aPairElementId])
                                {
                                    processedElements[anElementId] = true; 
                                    processedElements[aPairElementId] = true;

                                    Integer aNewNodeId = aMesh.nextNodeId();
                                    aMesh.addNode(aNewNodeId, aMidNode);

                                    anElement.setNodeId(0, aNodeId3);
                                    anElement.setNodeId(1, aNodeId1);
                                    anElement.setNodeId(2, aNewNodeId);

                                    aPairElement.setNodeId(0, aNodeId1);
                                    aPairElement.setNodeId(1, aNodeId4);
                                    aPairElement.setNodeId(2, aNewNodeId);

                                    aFace.setNodeId(0, aNodeId1);
                                    aFace.setNodeId(1, aNewNodeId);

                                    aPairFace.setNodeId(0, aNewNodeId);
                                    aPairFace.setNodeId(1, aNodeId1);

                                    // New elements
                                    CMshElement<Real> aNewElement1(ET_TRIANGLE);
                                    aNewElement1.addNodeId(aNodeId2);
                                    aNewElement1.addNodeId(aNodeId3);
                                    aNewElement1.addNodeId(aNewNodeId);

                                    Integer aNewElementId1 = aMesh.nextElementId();
                                    aMesh.addElement(aNewElementId1, aNewElement1);

                                    CMshElement<Real> aNewElement2(ET_TRIANGLE);
                                    aNewElement2.addNodeId(aNodeId4);
                                    aNewElement2.addNodeId(aNodeId2);
                                    aNewElement2.addNodeId(aNewNodeId);

                                    Integer aNewElementId2 = aMesh.nextElementId();
                                    aMesh.addElement(aNewElementId2, aNewElement2);
                                }
                            }
                        }
                    }
                }
            }

            Integer aFirstFaceId = 0;

            for (Integer i = 0; i < aMesh.nbFaces(); ++i)
            {
                Integer aFaceId = aMesh.faceId(i);

                if (i == 0)
                    aFirstFaceId = aFaceId;

                CMshFace<Real>& aFace = aMesh.face(aFaceId);

                if (aFace.pairFaceId() >= aMesh.nextFaceId())
                    aFace.setPairFaceId(aFirstFaceId);
            }

            aMesh.generateFaces(aTolerance);

            return true;
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::generate(const ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, Real meshSize, Real minMeshSize, Real maxMeshSize, const Real aTolerance)
        {
            std::vector<CGeoCoordinate<Real>> sInteriorPoints;

            ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
            aAnaFunction.set(meshSize);

            return this->generate(anEdgeMesh, maxNbElements, sInteriorPoints, aAnaFunction, minMeshSize, maxMeshSize, aTolerance);
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::generate(const CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, Real meshSize, Real minMeshSize, Real maxMeshSize, const Real aTolerance)
        {
            ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
            aAnaFunction.set(meshSize);

            return this->generate(anEdgeMesh, maxNbElements, sInteriorPoints, aAnaFunction, minMeshSize, maxMeshSize, aTolerance);
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::generate(const CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minMeshSize, Real maxMeshSize, const Real aTolerance)
        {
            m_previousNbElements = 0;

            m_bStop = false;

            // Add boundary nodes to surface mesh
            m_surfaceMesh.reset();
            m_boundingBox.reset();

            for (Integer i = 0; i < anEdgeMesh.nbNodes(); ++i)
            {
                Integer aNodeId = anEdgeMesh.nodeId(i);
                CMshNode<Real> aNode = anEdgeMesh.node(aNodeId);

                m_surfaceMesh.addNode(aNodeId, aNode);

                // Add to bounding box
                m_boundingBox.addCoordinate(aNode);
            }

            m_boundingBox.grow(aTolerance);

            // Add boundary to advancing front and rtree
            m_anAdvFront.clear();

            m_anAdvFront.reserve(anEdgeMesh.nbElements() * 50);

            m_tree.reset();

            std::map<Integer, Integer> newEdgeIds;

            for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i)
            {
                Integer anAdvEdgeId = anEdgeMesh.elementId(i);
                const CMshElement<Real>& anElement = anEdgeMesh.element(anAdvEdgeId);

                if (anElement.elementType() == ET_BEAM)
                {
                    newEdgeIds[anAdvEdgeId] = i;
                }
            }

            Real smallestEdgeLength = std::numeric_limits<Real>::max();
            Integer aFirstIndex = 0;

            m_nextEdgeId = 0;

            for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i)
            {
                Integer anAdvEdgeId = anEdgeMesh.elementId(i);

                const CMshElement<Real>& anElement = anEdgeMesh.element(anAdvEdgeId);

                if (anElement.elementType() == ET_BEAM)
                {
                    SMshAdvancingFrontEdge<Real> anAdvEdge;

                    anAdvEdge.id = m_nextEdgeId++;
                    anAdvEdge.remove = false;
                    anAdvEdge.boundary = true;

                    for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                        anAdvEdge.nodeId[j] = anElement.nodeId(j);

                    anAdvEdge.elementId = std::numeric_limits<Integer>::max();
                    anAdvEdge.nodeNotId3 = std::numeric_limits<Integer>::max();

                    anAdvEdge.build(anEdgeMesh);

                    for (Integer j = 0; j < anElement.nbFaceIds(); ++j)
                    {
                        anAdvEdge.neighborId[j] = std::numeric_limits<Integer>::max();

                        Integer aFaceId = anElement.faceId(j);
                        const CMshFace<Real>& aFace = anEdgeMesh.face(aFaceId);

                        if (aFace.hasPair())
                        {
                            Integer aPairFaceId = aFace.pairFaceId();
                            Integer anElementId = anEdgeMesh.face(aPairFaceId).elementId();

                            if (newEdgeIds.find(anElementId) != newEdgeIds.end())
                            {
                                anAdvEdge.neighborId[j] = newEdgeIds.at(anElementId);
                            }
                            else
                                std::cout << "Error: element id = " << anElementId << " not found!" << std::endl;
                        }
                        else
                        {
                            throw std::runtime_error("Boundary is open!");
                        }
                    }

                    if (anAdvEdge.line.length() < smallestEdgeLength)
                    {
                        smallestEdgeLength = anAdvEdge.line.length();
                        aFirstIndex = static_cast<Integer>(m_anAdvFront.size());
                    }

                    m_anAdvFront.push_back(anAdvEdge);

                    this->addEdgeToRtree(anAdvEdge, aTolerance);
                }
            }

            // Add interior nodes
            for (Integer i = 0; i < static_cast<Integer>(sInteriorPoints.size()); ++i)
            {
                Integer aNewNodeId = m_surfaceMesh.nextNodeId();
                CMshNode<Real> aNewNode = sInteriorPoints[i];

                m_surfaceMesh.addNode(aNewNodeId, aNewNode);

                SNode<Real> anInteriorNode;

                anInteriorNode.id = i;
                anInteriorNode.remove = false;

                anInteriorNode.nodeId = aNewNodeId;

                m_interiorNodes.push_back(anInteriorNode);
            }

            // Start meshing interior
            Integer maxElem = maxNbElements;

            this->advancingFrontTriMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 0.75, 0.05, false, false, aTolerance);
            this->advancingFrontTriMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 0.75, 0.02, true, false, aTolerance);
            this->advancingFrontTriMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 1.50, 0.00, true, false, aTolerance);

            if (this->frontSize() > 0)
            {
                std::cout << "Meshing error!" << std::endl;
                throw(m_anAdvFront);
            }

            m_surfaceMesh.removeDanglingNodes();
            m_surfaceMesh.renumber();
            m_surfaceMesh.generateFaces(aTolerance);

            return true;
        }

        template <typename Real>
        Integer CMshTriangleMesher<Real>::frontSize()
        {
            Integer n = 0;

            for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i)
            {
                if (!m_anAdvFront[i].remove)
                    n++;
            }

            return n;
        }

        template <typename Real>
        bool CMshTriangleMesher<Real>::advancingFrontTriMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real minMeshSize, Real maxMeshSize, Real sizeFactor, Real shrinkFactor, Real expandFactor, Real minQuality, const bool bAddNodes, const bool bCheckDelaunay, const Real aTolerance)
        {
            bool res = true;

            Real x, y;

            meshSizeFunc.removeAllVariables();

            meshSizeFunc.defineVariable("x", x);
            meshSizeFunc.defineVariable("y", y);

            static const Real pi = std::acos(-1.0);

            for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i)
            {
                this->checkUpdate();

                if (m_bStop)
                    return false;

                if (maxNbElements > 0)
                {
                    if (m_surfaceMesh.nbElements() >= maxNbElements)
                    {
                        std::cout << "Max number of elements (" << maxNbElements << ") reached!" << std::endl;
                        return false;
                    }
                }

                if (!m_anAdvFront[i].remove)
                {
                    Integer anAdvEdgeId = i;

                    SMshAdvancingFrontEdge<Real>& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                    Integer aNodeId1 = anAdvEdge.nodeId[0];
                    Integer aNodeId2 = anAdvEdge.nodeId[1];

                    CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
                    CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);

                    Integer aNodeId3 = m_anAdvFront[anAdvEdge.neighborId[0]].nodeId[0];
                    Integer aNodeId4 = m_anAdvFront[anAdvEdge.neighborId[1]].nodeId[1];

                    CMshNode<Real> aMidNode = (aNode1 + aNode2) * 0.5;

                    x = aMidNode.x();
                    y = aMidNode.y();

                    CGeoVector<Real> a = aNode2 - aNode1;

                    Real requiredMeshSize = meshSizeFunc.evaluate();
                    Real localMeshSize = static_cast<Real>(a.norm());

                    Real aFactor = std::max<Real>(std::min<Real>(requiredMeshSize / (localMeshSize + aTolerance * aTolerance), 2.0), 0.5);

                    localMeshSize *= aFactor;
                    localMeshSize = std::max(localMeshSize, minMeshSize);
                    localMeshSize = std::min(localMeshSize, maxMeshSize);

                    Real baseHeightSize = localMeshSize * sqrt(3.0) / 2.0; // Equilateral triangle (height to edge ratio)

                    // Rotate vector by 90 degrees
                    CGeoVector<Real> v = a;
                    v.rotate(pi * 0.5);
                    v.normalize();

                    // Add point to form triangle with correct spacing
                    CMshNode<Real> aNewNode = aMidNode + v * baseHeightSize * sizeFactor;
                    Integer aNewNodeId = m_surfaceMesh.nextNodeId();

                    // Get closest edges
                    CGeoBoundingBox<Real> aBoundingBox;
                    aBoundingBox.addCoordinate(aNode1);
                    aBoundingBox.addCoordinate(aNode2);
                    aBoundingBox.addCoordinate(aNewNode);
                    aBoundingBox.grow(baseHeightSize * sizeFactor * 0.5);

                    std::vector<Integer> sEdges;
                    m_tree.find(sEdges, aBoundingBox);

                    sEdges.erase(std::remove(sEdges.begin(), sEdges.end(), anAdvEdgeId), sEdges.end());

                    // Check if a node exists in proximity
                    std::vector<Integer> sNodes;
                    this->findClosestNodes(sEdges, sNodes);

                    // Meshing priority
                    // Priority = 1: close hole
                    // Priority = 2: other nodes in vicinity (3, 4, other node)
                    // Priority = 3: new node forming correct spacing

                    if (aNodeId3 == aNodeId4)
                    {
                        CMshNode<Real>& aNode3 = m_surfaceMesh.node(aNodeId3);

                        CMshTriangle<Real> aNewTriangle;

                        aNewTriangle.addVertex(aNode1);
                        aNewTriangle.addVertex(aNode2);
                        aNewTriangle.addVertex(aNode3);

                        aNewTriangle.calculateArea();

                        if (aNewTriangle.normal().z() > aTolerance)
                        {
                            Integer wNodeId;

                            if (this->edgeOk(anAdvEdge, aNode1, aNode3, sEdges, aTolerance) && this->edgeOk(anAdvEdge, aNode2, aNode3, sEdges, aTolerance) && !this->triangleContainsNode(aNode1, aNode2, aNode3, wNodeId, sNodes, aTolerance))
                            {
                                this->addTriangle(anAdvEdge, aNodeId3, sEdges, aTolerance);
                                res = true;
                                continue;
                            }
                        }
                    }

                    Real qmax = 0.0;

                    // If a node exists snap to that node
                    CMshNode<Real> anExistingNode;
                    Integer anExistingNodeId = std::numeric_limits<Integer>::max();

                    for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j)
                    {
                        Integer aNodeId = sNodes[j];

                        if (aNodeId == aNodeId1 || aNodeId == aNodeId2)
                            continue;

                        CMshNode<Real>& aNode = m_surfaceMesh.node(aNodeId);

                        Real d = (aNode - aNewNode).norm();

                        Real factor = 1.0;

                        if (aNodeId == aNodeId3)
                        {
                            Real angle = anAdvEdge.line.vector().angle(m_anAdvFront[anAdvEdge.neighborId[0]].line.vector());

                            if (angle > pi * 0.5)
                                factor = 1.5;
                        }

                        if (aNodeId == aNodeId4)
                        {
                            Real angle = anAdvEdge.line.vector().angle(m_anAdvFront[anAdvEdge.neighborId[1]].line.vector());

                            if (angle > pi * 0.5)
                                factor = 1.5;
                        }

                        // Use closest node
                        if (d < baseHeightSize * sizeFactor * expandFactor * factor)
                        {
                            CMshTriangle<Real> aNewTriangle1;

                            aNewTriangle1.addVertex(aNode1);
                            aNewTriangle1.addVertex(aNode2);
                            aNewTriangle1.addVertex(aNode);

                            aNewTriangle1.calculateArea();

                            if (aNewTriangle1.normal().z() > aTolerance)
                            {
                                aNewTriangle1.calculateQuality();

                                Real q1 = aNewTriangle1.quality();

                                if (q1 > qmax)
                                    q1 += this->edgeOk(anAdvEdge, aNode1, aNode, sEdges, aTolerance) ? 0.0 : -2.0;

                                if (q1 > qmax)
                                    q1 += this->edgeOk(anAdvEdge, aNode2, aNode, sEdges, aTolerance) ? 0.0 : -2.0;

                                if (q1 > qmax)
                                {
                                    Integer wNodeId;
                                    q1 += !this->triangleContainsNode(aNode1, aNode2, aNode, wNodeId, sNodes, aTolerance) ? 0.0 : -2.0;
                                }

                                if (q1 > qmax)
                                {
                                    qmax = q1;
                                    anExistingNodeId = aNodeId;
                                    anExistingNode = aNode;
                                }
                            }
                        }
                    }

                    if (qmax > minQuality)
                    {
                        for (Integer k = 0; k < static_cast<Integer>(m_interiorNodes.size()); ++k)
                        {
                            if (m_interiorNodes[k].remove)
                                continue;

                            if (m_interiorNodes[k].nodeId == anExistingNodeId)
                                m_interiorNodes[k].remove = true;
                        }

                        this->addTriangle(anAdvEdge, anExistingNodeId, sEdges, aTolerance);
                        res = true;
                    }
                    else if (bAddNodes)
                    {
                        CGeoLine<Real> aLine(aMidNode, aNewNode);

                        Real dmin = findShortestDistance(sEdges, aLine, anAdvEdgeId, aTolerance);

                        if (dmin > baseHeightSize * sizeFactor * shrinkFactor * 0.25 && dmin < baseHeightSize * sizeFactor * expandFactor)
                            aNewNode = aMidNode + v * dmin * 0.5;

                        CMshTriangle<Real> aNewTriangle2;

                        aNewTriangle2.addVertex(aNode1);
                        aNewTriangle2.addVertex(aNode2);
                        aNewTriangle2.addVertex(aNewNode);

                        aNewTriangle2.calculateArea();

                        if (aNewTriangle2.normal().z() > aTolerance)
                        {
                            aNewTriangle2.calculateQuality();

                            Real q2 = aNewTriangle2.quality();

                            if (q2 > 1.5 * minQuality)
                                q2 += this->edgeOk(anAdvEdge, aNode1, aNewNode, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q2 > 1.5 * minQuality)
                                q2 += this->edgeOk(anAdvEdge, aNode2, aNewNode, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q2 > 1.5 * minQuality)
                            {
                                Integer wNodeId;
                                q2 += !this->triangleContainsNode(aNode1, aNode2, aNewNode, wNodeId, sNodes, aTolerance) ? 0.0 : -2.0;
                            }

                            if (q2 > 1.5 * minQuality && bCheckDelaunay)
                                q2 += this->checkDelaunay(aNewNode, aTolerance) ? 0.0 : -2.0;

                            if (q2 > 1.5 * minQuality)
                            {
                                // Create a new node
                                m_surfaceMesh.addNode(aNewNodeId, aNewNode);

                                this->addTriangle(anAdvEdge, aNewNodeId, sEdges, aTolerance);
                                res = true;

                                if (!m_boundingBox.contains(aNewNode, aTolerance))
                                {
                                    throw std::runtime_error("Node is outside boundary!");
                                }
                            }
                        }
                    }
                }
            }

            if (onUpdate != nullptr)
                onUpdate(0);

            return res;
        }

        template <typename Real>
        CMshMesh<Real>& CMshTriangleMesher<Real>::mesh()
        {
            return m_surfaceMesh;
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::setIntervals(const Integer timeInterval, const Integer dataInterval)
        {
            m_begin = clock();

            m_timeInterval = timeInterval;
            m_dataInterval = dataInterval;
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::stopMeshing()
        {
            m_bStop = true;
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::applyFixedBoundary(ENigMA::mesh::CMshMesh<Real>& aSurfaceMesh, CMshMesh<Real>& anEdgeMesh, const Real aTolerance)
        {
            // Discover double faces
            CGeoHashGrid<Real> aHashGrid;

            std::vector<CGeoCoordinate<Real>> sCenterCoordinates;

            for (Integer i = 0; i < aSurfaceMesh.nbFaces(); ++i)
            {
                Integer aFaceId = aSurfaceMesh.faceId(i);
                CMshFace<Real>& aFace = aSurfaceMesh.face(aFaceId);

                CGeoCoordinate<Real> aCenterCoordinate(0.0, 0.0, 0.0);

                for (Integer j = 0; j < aFace.nbNodeIds(); ++j)
                {
                    Integer aNodeId = aFace.nodeId(j);
                    CMshNode<Real>& aNode = aSurfaceMesh.node(aNodeId);

                    aCenterCoordinate += aNode;
                }

                if (aFace.nbNodeIds() > 0)
                    aCenterCoordinate /= static_cast<Real>(aFace.nbNodeIds());

                sCenterCoordinates.push_back(aCenterCoordinate);

                aHashGrid.addGeometricObject(aFaceId, aCenterCoordinate);
            }

            aHashGrid.build();

            for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i)
            {
                Integer anElementId = anEdgeMesh.elementId(i);

                CMshElement<Real>& anElement = anEdgeMesh.element(anElementId);

                CGeoCoordinate<Real> aCenterCoordinate(0.0, 0.0, 0.0);

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    Integer aNodeId = anElement.nodeId(j);
                    CMshNode<Real> aNode = anEdgeMesh.node(aNodeId);

                    aCenterCoordinate += aNode;
                }

                if (anElement.nbNodeIds() > 0)
                    aCenterCoordinate /= static_cast<Real>(anElement.nbNodeIds());

                std::vector<Integer> sCoordinates;

                aHashGrid.find(sCoordinates, aCenterCoordinate, aTolerance);

                for (Integer j = 0; j < static_cast<Integer>(sCoordinates.size()); ++j)
                {
                    Integer aFaceId = sCoordinates[j];
                    aSurfaceMesh.face(aFaceId).setHasPair(false);
                }
            }
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::flipEdges(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aTolerance)
        {
            std::map<Integer, bool> sFlipped;

            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {
                Integer anElementId = aMesh.elementId(i);
                sFlipped[anElementId] = false;
            }

            for (Integer i = 0; i < aMesh.nbFaces(); ++i)
            {
                Integer aFaceId = aMesh.faceId(i);
                CMshFace<Real> aFace = aMesh.face(aFaceId);

                if (aFace.faceType() != FT_LINE)
                    continue;

                Integer anElementId = aFace.elementId();

                if (sFlipped.at(anElementId))
                    continue;

                CMshElement<Real>& anElement = aMesh.element(anElementId);

                if (anElement.elementType() != ET_TRIANGLE)
                    continue;

                Integer aNodeId1 = aFace.nodeId(0);
                Integer aNodeId2 = aFace.nodeId(1);
                Integer aNodeId3;
                Integer aNodeId4;

                bool bNodeId3 = false;
                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    if (anElement.nodeId(j) != aNodeId1 && anElement.nodeId(j) != aNodeId2)
                    {
                        aNodeId3 = anElement.nodeId(j);
                        bNodeId3 = true;
                        break;
                    }
                }

                if (!bNodeId3)
                    continue;

                if (aFace.hasPair())
                {
                    Integer pairFaceId = aFace.pairFaceId();
                    CMshFace<Real> aPairFace = aMesh.face(pairFaceId);

                    Integer aNeighborId = aPairFace.elementId();

                    if (sFlipped.at(aNeighborId))
                        continue;

                    CMshElement<Real>& aNeighbor = aMesh.element(aNeighborId);

                    if (aNeighbor.elementType() != ET_TRIANGLE)
                        continue;

                    bool bNodeId4 = false;
                    for (Integer j = 0; j < aNeighbor.nbNodeIds(); ++j)
                    {
                        if (aNeighbor.nodeId(j) != aNodeId1 && aNeighbor.nodeId(j) != aNodeId2)
                        {
                            aNodeId4 = aNeighbor.nodeId(j);
                            bNodeId4 = true;
                            break;
                        }
                    }

                    if (!bNodeId4)
                        continue;

                    const CMshNode<Real>& aNode1 = aMesh.node(aNodeId1);
                    const CMshNode<Real>& aNode2 = aMesh.node(aNodeId2);
                    const CMshNode<Real>& aNode3 = aMesh.node(aNodeId3);
                    const CMshNode<Real>& aNode4 = aMesh.node(aNodeId4);

                    // Original
                    CMshTriangle<Real> aTriangle1;
                    aTriangle1.addVertex(aNode1);
                    aTriangle1.addVertex(aNode2);
                    aTriangle1.addVertex(aNode3);

                    aTriangle1.calculateQuality();
                    Real q1 = aTriangle1.quality();

                    CMshTriangle<Real> aTriangle2;
                    aTriangle2.addVertex(aNode2);
                    aTriangle2.addVertex(aNode1);
                    aTriangle2.addVertex(aNode4);

                    aTriangle2.calculateQuality();
                    Real q2 = aTriangle2.quality();

                    // Modified
                    CMshTriangle<Real> aTriangle3;
                    aTriangle3.addVertex(aNode3);
                    aTriangle3.addVertex(aNode4);
                    aTriangle3.addVertex(aNode2);

                    aTriangle3.calculateArea();
                    aTriangle3.calculateQuality();
                    Real q3 = aTriangle3.quality();

                    CMshTriangle<Real> aTriangle4;
                    aTriangle4.addVertex(aNode4);
                    aTriangle4.addVertex(aNode3);
                    aTriangle4.addVertex(aNode1);

                    aTriangle4.calculateArea();
                    aTriangle4.calculateQuality();
                    Real q4 = aTriangle4.quality();

                    if (std::min(q3, q4) > std::min(q1, q2) && aTriangle3.normal().z() > aTolerance && aTriangle4.normal().z() > aTolerance)
                    {
                        // Do flip
                        aMesh.element(anElementId).setNodeId(0, aNodeId3);
                        aMesh.element(anElementId).setNodeId(1, aNodeId4);
                        aMesh.element(anElementId).setNodeId(2, aNodeId2);

                        aMesh.element(aNeighborId).setNodeId(0, aNodeId4);
                        aMesh.element(aNeighborId).setNodeId(1, aNodeId3);
                        aMesh.element(aNeighborId).setNodeId(2, aNodeId1);

                        sFlipped[anElementId] = true;
                        sFlipped[aNeighborId] = true;
                    }
                }
            }

            aMesh.generateFaces(aTolerance);
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::relaxNodes(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aTolerance)
        {
            std::map<Integer, bool> bBoundaryNode;

            for (Integer i = 0; i < aMesh.nbNodes(); ++i)
            {
                Integer aNodeId = aMesh.nodeId(i);
                bBoundaryNode[aNodeId] = false;
            }

            for (Integer i = 0; i < aMesh.nbFaces(); ++i)
            {
                Integer aFaceId = aMesh.faceId(i);
                CMshFace<Real>& aFace = aMesh.face(aFaceId);

                if (aFace.hasPair())
                    continue;

                for (Integer j = 0; j < aFace.nbNodeIds(); ++j)
                {
                    Integer aNodeId = aFace.nodeId(j);
                    bBoundaryNode[aNodeId] = true;
                }
            }

            CMshMeshQuery<Real> aMeshQuery(aMesh);

            std::vector<Integer> sElementIds;

            for (Integer i = 0; i < aMesh.nbNodes(); ++i)
            {
                Integer aMovingNodeId = aMesh.nodeId(i);
                CMshNode<Real>& aMovingNode = aMesh.node(aMovingNodeId);

                if (bBoundaryNode.at(aMovingNodeId))
                    continue;

                aMeshQuery.elementsSharingNode(aMovingNodeId, sElementIds);

                std::vector<Integer> sNodeIds;
                sNodeIds.push_back(aMovingNodeId);

                for (Integer j = 0; j < static_cast<Integer>(sElementIds.size()); ++j)
                {
                    Integer anElementId = sElementIds[j];

                    CMshElement<Real>& anElement = aMesh.element(anElementId);

                    if (anElement.elementType() != ET_TRIANGLE)
                        continue;

                    for (Integer k = 0; k < anElement.nbNodeIds(); ++k)
                    {
                        Integer aNodeId = anElement.nodeId(k);

                        if (std::find(sNodeIds.begin(), sNodeIds.end(), aNodeId) == sNodeIds.end())
                            sNodeIds.push_back(aNodeId);
                    }
                }

                if (sNodeIds.size() > 1)
                {
                    CMshNode<Real> aNewNode(0.0, 0.0, 0.0);

                    for (Integer k = 0; k < static_cast<Integer>(sNodeIds.size()); ++k)
                        aNewNode += aMesh.node(sNodeIds[k]);

                    aNewNode /= static_cast<Real>(sNodeIds.size());

                    Real sumq1 = 0.0;
                    Real sumq2 = 0.0;

                    bool bMoveNode = true;

                    for (Integer k = 0; k < static_cast<Integer>(sElementIds.size()); ++k)
                    {
                        CMshElement<Real>& aModElement = aMesh.element(sElementIds[k]);

                        if (aModElement.elementType() != ET_TRIANGLE)
                            continue;

                        Integer aNodeId1 = aModElement.nodeId(0);
                        Integer aNodeId2 = aModElement.nodeId(1);
                        Integer aNodeId3 = aModElement.nodeId(2);

                        CMshNode<Real>& aNode1 = aMesh.node(aNodeId1);
                        CMshNode<Real>& aNode2 = aMesh.node(aNodeId2);
                        CMshNode<Real>& aNode3 = aMesh.node(aNodeId3);

                        CGeoCoordinate<Real> aVertex1 = aNode1;
                        CGeoCoordinate<Real> aVertex2 = aNode2;
                        CGeoCoordinate<Real> aVertex3 = aNode3;

                        CMshTriangle<Real> aTriangle1;

                        aTriangle1.addVertex(aVertex1);
                        aTriangle1.addVertex(aVertex2);
                        aTriangle1.addVertex(aVertex3);

                        aTriangle1.calculateQuality();

                        Real q1 = aTriangle1.quality();

                        if (aNodeId1 == aMovingNodeId)
                            aVertex1 = aNewNode;

                        if (aNodeId2 == aMovingNodeId)
                            aVertex2 = aNewNode;

                        if (aNodeId3 == aMovingNodeId)
                            aVertex3 = aNewNode;

                        CMshTriangle<Real> aTriangle2;

                        aTriangle2.addVertex(aVertex1);
                        aTriangle2.addVertex(aVertex2);
                        aTriangle2.addVertex(aVertex3);

                        aTriangle2.calculateArea();
                        aTriangle2.calculateQuality();

                        Real q2 = aTriangle2.quality();

                        sumq1 += q1;
                        sumq2 += q2;

                        if (aTriangle2.normal().z() <= aTolerance)
                        {
                            bMoveNode = false;
                            break;
                        }
                    }

                    if (bMoveNode && sumq2 > sumq1)
                        aMovingNode = aNewNode;
                }
            }
        }

        template <typename Real>
        void CMshTriangleMesher<Real>::collapseEdges(ENigMA::mesh::CMshMesh<Real>& aMesh, Real collapseSize, const Real aTolerance)
        {
            std::map<Integer, bool> bBoundaryNode;

            for (Integer i = 0; i < aMesh.nbNodes(); ++i)
            {
                Integer aNodeId = aMesh.nodeId(i);
                bBoundaryNode[aNodeId] = false;
            }

            for (Integer i = 0; i < aMesh.nbFaces(); ++i)
            {
                Integer aFaceId = aMesh.faceId(i);
                CMshFace<Real>& aFace = aMesh.face(aFaceId);

                if (aFace.hasPair())
                    continue;

                for (Integer j = 0; j < aFace.nbNodeIds(); ++j)
                {
                    Integer aNodeId = aFace.nodeId(j);
                    bBoundaryNode[aNodeId] = true;
                }
            }

            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {
                Integer anElementId = aMesh.elementId(i);
                CMshElement<Real>& anElement = aMesh.element(anElementId);

                if (anElement.elementType() != ET_TRIANGLE)
                    continue;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    Integer aNodeId1 = anElement.nodeId((j + 0) % 3);
                    Integer aNodeId2 = anElement.nodeId((j + 1) % 3);

                    CMshNode<Real>& aNode1 = aMesh.node(aNodeId1);
                    CMshNode<Real>& aNode2 = aMesh.node(aNodeId2);

                    Real anEdgeLength = (aNode1 - aNode2).norm();

                    if (bBoundaryNode.at(aNodeId2))
                        continue;

                    if (anEdgeLength < collapseSize * 0.75)
                    {
                        // Try to collapse
                        aNode2 = aNode1;
                        anElement.setElementType(ET_NONE);
                        break;
                    }
                }
            }

            // Delete invalid elements
            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {
                Integer anElementId = aMesh.elementId(i);
                CMshElement<Real>& anElement = aMesh.element(anElementId);

                if (anElement.elementType() == ET_NONE)
                {
                    aMesh.removeElement(anElementId);
                }
            }

            aMesh.mergeNodes(aTolerance);
            aMesh.renumber();

            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {
                Integer anElementId = aMesh.elementId(i);
                CMshElement<Real>& anElement = aMesh.element(anElementId);

                if (anElement.elementType() == ET_NONE)
                {
                    aMesh.element(anElementId).setElementType(ET_TRIANGLE);
                }
            }

            aMesh.generateFaces(aTolerance);
        }
    }
}
