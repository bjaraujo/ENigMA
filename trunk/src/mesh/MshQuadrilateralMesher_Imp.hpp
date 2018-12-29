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

#include "GeoLine.hpp"
#include "GeoCircle.hpp"
#include "MshTriangle.hpp"
#include "MshQuadrilateral.hpp"

namespace ENigMA
{

    namespace mesh
    {

        template <typename Real>
        CMshQuadrilateralMesher<Real>::CMshQuadrilateralMesher() : 
            m_begin(), m_end(), m_nextEdgeId(0),
            m_previousNbElements(0), 
            m_timeInterval(1), m_dataInterval(100), m_bStop(false), onUpdate(nullptr)
        {

        }

        template <typename Real>
        CMshQuadrilateralMesher<Real>::~CMshQuadrilateralMesher()
        {

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::checkUpdate()
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
        void CMshQuadrilateralMesher<Real>::removeEdge(SAdvancingFrontEdge& anAdvEdge, const Real aTolerance)
        {

            anAdvEdge.remove = true;
            removeEdgeFromRtree(anAdvEdge, aTolerance);

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::addEdgeToRtree(SAdvancingFrontEdge& anAdvEdge, const Real aTolerance)
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
        void CMshQuadrilateralMesher<Real>::removeEdgeFromRtree(SAdvancingFrontEdge& anAdvEdge, const Real aTolerance)
        {

            Integer aNodeId1 = anAdvEdge.nodeId[0];
            Integer aNodeId2 = anAdvEdge.nodeId[1];

            CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
            CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);

            CGeoBoundingBox<Real> aBoundingBox;

            aBoundingBox.addCoordinate(aNode1);
            aBoundingBox.addCoordinate(aNode2);

            aBoundingBox.grow(aTolerance);

            m_tree.removeGeometricObject(anAdvEdge.id, aBoundingBox);

        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::edgeExists(SAdvancingFrontEdge& anAdvEdge, Integer& aDuplicateEdgeId, std::vector<Integer>& sEdges, const Real aTolerance)
        {

            std::vector<Integer> sNodeIds;

            sNodeIds.push_back(anAdvEdge.nodeId[0]);
            sNodeIds.push_back(anAdvEdge.nodeId[1]);

            std::sort(sNodeIds.begin(), sNodeIds.end());

            for (Integer j = 0; j < static_cast<Integer> (sEdges.size()); ++j)
            {

                // Discard same edge
                if (sEdges[j] == anAdvEdge.id)
                    continue;

                SAdvancingFrontEdge& anotherEdge = m_anAdvFront[sEdges[j]];

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
        bool CMshQuadrilateralMesher<Real>::edgeOk(SAdvancingFrontEdge& anAdvEdge, CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, std::vector<Integer>& sEdges, const Real aTolerance)
        {

            CGeoCoordinate<Real> aPoint;

            CGeoLine<Real> aLine1(aNode1, aNode2);

            for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j)
            {

                SAdvancingFrontEdge& anotherEdge = m_anAdvFront[sEdges[j]];

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
        void CMshQuadrilateralMesher<Real>::findClosestNodes(std::vector<Integer>& sEdges, std::vector<Integer>& sNodes)
        {

            sNodes.clear();

            for (Integer j = 0; j < static_cast<Integer> (sEdges.size()); ++j)
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

            for (Integer j = 0; j < static_cast<Integer> (m_interiorNodes.size()); ++j)
            {

                if (m_interiorNodes[j].remove)
                    continue;

                Integer aNodeId = m_interiorNodes[j].nodeId;

                if (std::find(sNodes.begin(), sNodes.end(), aNodeId) == sNodes.end())
                    sNodes.push_back(aNodeId);

            }

        }

        template <typename Real>
        Real CMshQuadrilateralMesher<Real>::findShortestDistance(std::vector<Integer>& sEdges, CGeoLine<Real>& aLine, Integer anAdvEdgeId, const Real aTolerance)
        {

            Real dmin = std::numeric_limits<Real>::max();

            for (Integer j = 0; j < static_cast<Integer> (sEdges.size()); ++j)
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
        bool CMshQuadrilateralMesher<Real>::quadrilateralContainsNode(CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, CMshNode<Real>& aNode4, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance)
        {

            CMshQuadrilateral<Real> aQuadrilateral;

            aQuadrilateral.addVertex(aNode1);
            aQuadrilateral.addVertex(aNode2);
            aQuadrilateral.addVertex(aNode3);
            aQuadrilateral.addVertex(aNode4);

            for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j)
            {

                aNodeId = sNodes[j];

                CMshNode<Real>& aNode = m_surfaceMesh.node(aNodeId);

                if ((aNode - aNode1).norm() < aTolerance ||
                    (aNode - aNode2).norm() < aTolerance ||
                    (aNode - aNode3).norm() < aTolerance ||
                    (aNode - aNode4).norm() < aTolerance)
                    continue;

                CGeoIntersectionType anIntersectionType;

                if (aQuadrilateral.contains(aNode, anIntersectionType, aTolerance))
                {

                    if (anIntersectionType == IT_INTERNAL || anIntersectionType == IT_EDGE)
                        return true;

                }

            }

            return false;

        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::checkDelaunay(CMshNode<Real>& aNewNode, const Real aTolerance)
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
        void CMshQuadrilateralMesher<Real>::adjustConnectivity(std::vector<Integer>& sEdges)
        {

            try
            {

                for (Integer i = 0; i < static_cast<Integer> (m_anAdvFront.size()); ++i)
                {

                    Integer anAdvEdgeId = i;

                    SAdvancingFrontEdge& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                    if (anAdvEdge.remove)
                        continue;

                    if (anAdvEdge.nodeId[0] == anAdvEdge.nodeId[1])
                    {
                        anAdvEdge.remove = true;
                        continue;
                    }

                    if (anAdvEdge.neighborId[0] >= static_cast<Integer>(m_anAdvFront.size()) ||
                        anAdvEdge.neighborId[1] >= static_cast<Integer>(m_anAdvFront.size()))
                    {
                        throw std::out_of_range("Connectivity is out of range!");
                    }

                    SAdvancingFrontEdge& aPrevEdge = m_anAdvFront[anAdvEdge.neighborId[0]];
                    SAdvancingFrontEdge& aNextEdge = m_anAdvFront[anAdvEdge.neighborId[1]];

                    if (aPrevEdge.remove || aNextEdge.remove)
                    {

                        for (Integer j = 0; j < static_cast<Integer> (sEdges.size()); ++j)
                        {

                            Integer anotherEdgeId = sEdges[j];

                            // Discard same edge
                            if (anotherEdgeId == anAdvEdge.id)
                                continue;

                            SAdvancingFrontEdge& anotherEdge = m_anAdvFront[anotherEdgeId];

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
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::cleanDuplicateEdges(std::vector<Integer>& sEdges, const Real aTolerance)
        {

            try
            {

                if (m_nextEdgeId < 3)
                    return;

                bool bAdjustConnectivity = false;

                // Check last three edges
                for (Integer i = 0; i < 3; ++i)
                {

                    Integer anAdvEdgeId = m_nextEdgeId - i - 1;

                    SAdvancingFrontEdge& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                    if (anAdvEdge.remove)
                        continue;

                    Integer aDuplicateEdgeId;

                    if (edgeExists(anAdvEdge, aDuplicateEdgeId, sEdges, aTolerance))
                    {

                        SAdvancingFrontEdge& aDuplicateEdge = m_anAdvFront[aDuplicateEdgeId];

                        this->removeEdge(anAdvEdge, aTolerance);
                        this->removeEdge(aDuplicateEdge, aTolerance);

                        bAdjustConnectivity = true;

                    }

                }

                if (bAdjustConnectivity)
                    this->adjustConnectivity(sEdges);

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::addQuadrilateral(SAdvancingFrontEdge& anAdvEdge, const Integer aNodeId3, const Integer aNodeId4, std::vector<Integer>& sEdges, const Real aTolerance)
        {

            try
            {

                Integer aNodeId1 = anAdvEdge.nodeId[0];
                Integer aNodeId2 = anAdvEdge.nodeId[1];

                if (aNodeId3 == aNodeId1 ||
                    aNodeId3 == aNodeId2 ||
                    aNodeId4 == aNodeId1 ||
                    aNodeId4 == aNodeId2)
                {
                    std::cout << "Error: invalid quadrilateral!" << std::endl;
                    return;
                }

                if (aNodeId3 == aNodeId4)
                {

                    // Add new triangle
                    CMshElement<Real> aNewElement(ET_TRIANGLE);
                    aNewElement.addNodeId(aNodeId1);
                    aNewElement.addNodeId(aNodeId2);
                    aNewElement.addNodeId(aNodeId3);

                    Integer aNewElementId = m_surfaceMesh.nextElementId();

                    m_surfaceMesh.addElement(aNewElementId, aNewElement);

                    SAdvancingFrontEdge& aPrevEdge = m_anAdvFront[anAdvEdge.neighborId[0]];
                    SAdvancingFrontEdge& aNextEdge = m_anAdvFront[anAdvEdge.neighborId[1]];

                    Integer aNewEdgeId1 = m_nextEdgeId++;
                    Integer aNewEdgeId2 = m_nextEdgeId++;

                    // Add edge 1
                    SAdvancingFrontEdge aNewEdge1;
                    aNewEdge1.id = aNewEdgeId1;
                    aNewEdge1.remove = false;
                    aNewEdge1.boundary = false;
                    aNewEdge1.nodeId[0] = aNodeId1;
                    aNewEdge1.nodeId[1] = aNodeId3;
                    aNewEdge1.neighborId[0] = anAdvEdge.neighborId[0];
                    aNewEdge1.neighborId[1] = aNewEdgeId2;
                    aNewEdge1.quadrilateralId = aNewElementId;

                    aNewEdge1.build(m_surfaceMesh);

                    // Correct connectivity
                    aPrevEdge.neighborId[1] = aNewEdgeId1;

                    // Add edge to rtree
                    addEdgeToRtree(aNewEdge1, aTolerance);

                    // Add edge 2
                    SAdvancingFrontEdge aNewEdge2;
                    aNewEdge2.id = aNewEdgeId2;
                    aNewEdge2.remove = false;
                    aNewEdge2.boundary = false;
                    aNewEdge2.nodeId[0] = aNodeId3;
                    aNewEdge2.nodeId[1] = aNodeId2;
                    aNewEdge2.neighborId[0] = aNewEdgeId1;
                    aNewEdge2.neighborId[1] = anAdvEdge.neighborId[1];
                    aNewEdge2.quadrilateralId = aNewElementId;

                    aNewEdge2.build(m_surfaceMesh);

                    // Correct connectivity
                    aNextEdge.neighborId[0] = aNewEdgeId2;

                    // Add edge to rtree
                    addEdgeToRtree(aNewEdge2, aTolerance);

                    // Remove this edge
                    this->removeEdge(anAdvEdge, aTolerance);

                    m_anAdvFront.push_back(aNewEdge1);
                    m_anAdvFront.push_back(aNewEdge2);

                    this->cleanDuplicateEdges(sEdges, aTolerance);

                }
                else
                {

                    // Add new quadrilateral
                    CMshElement<Real> aNewElement(ET_QUADRILATERAL);
                    aNewElement.addNodeId(aNodeId1);
                    aNewElement.addNodeId(aNodeId2);
                    aNewElement.addNodeId(aNodeId4);
                    aNewElement.addNodeId(aNodeId3);

                    Integer aNewElementId = m_surfaceMesh.nextElementId();

                    m_surfaceMesh.addElement(aNewElementId, aNewElement);

                    SAdvancingFrontEdge& aPrevEdge = m_anAdvFront[anAdvEdge.neighborId[0]];
                    SAdvancingFrontEdge& aNextEdge = m_anAdvFront[anAdvEdge.neighborId[1]];

                    Integer aNewEdgeId1 = m_nextEdgeId++;
                    Integer aNewEdgeId2 = m_nextEdgeId++;
                    Integer aNewEdgeId3 = m_nextEdgeId++;

                    // Add edge 1
                    SAdvancingFrontEdge aNewEdge1;
                    aNewEdge1.id = aNewEdgeId1;
                    aNewEdge1.remove = false;
                    aNewEdge1.boundary = false;
                    aNewEdge1.nodeId[0] = aNodeId1;
                    aNewEdge1.nodeId[1] = aNodeId3;
                    aNewEdge1.neighborId[0] = anAdvEdge.neighborId[0];
                    aNewEdge1.neighborId[1] = aNewEdgeId2;
                    aNewEdge1.quadrilateralId = aNewElementId;

                    aNewEdge1.build(m_surfaceMesh);

                    // Correct connectivity
                    aPrevEdge.neighborId[1] = aNewEdgeId1;

                    // Add edge to rtree
                    addEdgeToRtree(aNewEdge1, aTolerance);

                    // Add edge 2
                    SAdvancingFrontEdge aNewEdge2;
                    aNewEdge2.id = aNewEdgeId2;
                    aNewEdge2.remove = false;
                    aNewEdge2.boundary = false;
                    aNewEdge2.nodeId[0] = aNodeId3;
                    aNewEdge2.nodeId[1] = aNodeId4;
                    aNewEdge2.neighborId[0] = aNewEdgeId1;
                    aNewEdge2.neighborId[1] = aNewEdgeId3;
                    aNewEdge2.quadrilateralId = aNewElementId;

                    aNewEdge2.build(m_surfaceMesh);

                    // Add edge to rtree
                    addEdgeToRtree(aNewEdge2, aTolerance);

                    // Add edge 3
                    SAdvancingFrontEdge aNewEdge3;
                    aNewEdge3.id = aNewEdgeId3;
                    aNewEdge3.remove = false;
                    aNewEdge3.boundary = false;
                    aNewEdge3.nodeId[0] = aNodeId4;
                    aNewEdge3.nodeId[1] = aNodeId2;
                    aNewEdge3.neighborId[0] = aNewEdgeId2;
                    aNewEdge3.neighborId[1] = anAdvEdge.neighborId[1];
                    aNewEdge3.quadrilateralId = aNewElementId;

                    aNewEdge3.build(m_surfaceMesh);

                    // Correct connectivity
                    aNextEdge.neighborId[0] = aNewEdgeId3;

                    // Add edge to rtree
                    addEdgeToRtree(aNewEdge3, aTolerance);

                    // Remove this edge
                    this->removeEdge(anAdvEdge, aTolerance);

                    m_anAdvFront.push_back(aNewEdge1);
                    m_anAdvFront.push_back(aNewEdge2);
                    m_anAdvFront.push_back(aNewEdge3);

                    this->cleanDuplicateEdges(sEdges, aTolerance);

                }

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }

        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::remesh(CMshMesh<Real>& anEdgeMesh, Real meshSize)
        {

            std::stringstream ss(std::stringstream::in | std::stringstream::out);

            ss << meshSize;

            ENigMA::analytical::CAnaFunction<Real> aAnaFunction;

            aAnaFunction.set(ss.str());

            return this->remesh(anEdgeMesh, aAnaFunction);

        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::remesh(CMshMesh<Real>& anEdgeMesh, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc)
        {

            Real x, y;

            meshSizeFunc.removeAllVariables();

            meshSizeFunc.defineVariable("x", x);
            meshSizeFunc.defineVariable("y", y);

            // Split edge mesh according to local mesh size
            Integer nbEdges = anEdgeMesh.nbElements();

            for (Integer i = 0; i < nbEdges; ++i)
            {

                Integer anAdvEdgeId = anEdgeMesh.elementId(i);

                CMshElement<Real>& anElement = anEdgeMesh.element(anAdvEdgeId);

                Integer aNodeId1 = anElement.nodeId(0);
                Integer aNodeId2 = anElement.nodeId(1);

                CMshNode<Real> aNode1 = anEdgeMesh.node(aNodeId1);
                CMshNode<Real> aNode2 = anEdgeMesh.node(aNodeId2);

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

                        Integer aNewNodeId = anEdgeMesh.nextNodeId();
                        anEdgeMesh.addNode(aNewNodeId, aNode);

                        CMshElement<Real> aNewElement(ET_BEAM);
                        aNewElement.addNodeId(aPrevNodeId);
                        aNewElement.addNodeId(aNewNodeId);

                        Integer aNewElementId = anEdgeMesh.nextElementId();

                        if (bConnectivity)
                        {

                            std::vector<CMshFace<Real> > sFaces;
                            aNewElement.generateFaces(sFaces);

                            for (Integer j = 0; j < static_cast<Integer> (sFaces.size()); ++j)
                            {

                                CMshFace<Real> aNewFace = sFaces[j];

                                Integer aNewFaceId = anEdgeMesh.nextFaceId();

                                aNewFace.setElementId(aNewElementId);

                                if (j == 0)
                                    aNewFace.setPairFaceId(aPrevFaceId);
                                else
                                    aNewFace.setPairFaceId(aNewFaceId + 1);

                                anEdgeMesh.addFace(aNewFaceId, aNewFace);
                                aNewElement.addFaceId(aNewFaceId);

                                aPrevFaceId = aNewFaceId;

                            }

                        }

                        aPrevNodeId = aNewNodeId;

                        anEdgeMesh.addElement(aNewElementId, aNewElement);

                        // Split element
                        anEdgeMesh.element(anAdvEdgeId).setNodeId(0, aNewNodeId);
                        anEdgeMesh.element(anAdvEdgeId).setNodeId(1, aNodeId2);

                        if (bConnectivity)
                        {
                            anEdgeMesh.element(anAdvEdgeId).setFaceId(0, aPrevFaceId);
                            anEdgeMesh.element(anAdvEdgeId).setFaceId(1, aFaceId2);
                        }

                    }

                }

            }

            Integer aFirstFaceId = 0;

            for (Integer i = 0; i < anEdgeMesh.nbFaces(); ++i)
            {

                Integer aFaceId = anEdgeMesh.faceId(i);

                if (i == 0)
                    aFirstFaceId = aFaceId;

                CMshFace<Real>& aFace = anEdgeMesh.face(aFaceId);

                if (aFace.pairFaceId() >= anEdgeMesh.nextFaceId())
                    aFace.setPairFaceId(aFirstFaceId);

            }

            return true;

        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::generate(CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, Real meshSize, Real minQuality, const Real aTolerance)
        {

            try
            {

                std::vector<CGeoCoordinate<Real> > sInteriorPoints;

                ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
                aAnaFunction.set(meshSize);

                return this->generate(anEdgeMesh, maxNbElements, sInteriorPoints, aAnaFunction, minQuality, aTolerance);

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                return false;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                return false;
            }

        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::generate(CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real> >& sInteriorPoints, Real meshSize, Real minQuality, const Real aTolerance)
        {

            try
            {

                ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
                aAnaFunction.set(meshSize);

                return this->generate(anEdgeMesh, maxNbElements, sInteriorPoints, aAnaFunction, minQuality, aTolerance);

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                return false;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                return false;
            }

        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::generate(CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real> >& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minQuality, const Real aTolerance)
        {

            try
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

                m_tree.reset();

                m_nextEdgeId = 0;

                std::map<Integer, Integer> newEdgeIds;

                for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i)
                {

                    Integer anAdvEdgeId = anEdgeMesh.elementId(i);

                    newEdgeIds[anAdvEdgeId] = i;
                }

                for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i)
                {

                    Integer anAdvEdgeId = anEdgeMesh.elementId(i);

                    CMshElement<Real>& anElement = anEdgeMesh.element(anAdvEdgeId);

                    if (anElement.elementType() == ET_BEAM)
                    {

                        SAdvancingFrontEdge anAdvEdge;

                        anAdvEdge.id = m_nextEdgeId++;
                        anAdvEdge.remove = false;
                        anAdvEdge.boundary = true;

                        for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                            anAdvEdge.nodeId[j] = anElement.nodeId(j);

                        anAdvEdge.quadrilateralId = std::numeric_limits<Integer>::max();

                        anAdvEdge.build(anEdgeMesh);

                        for (Integer j = 0; j < anElement.nbFaceIds(); ++j)
                        {

                            anAdvEdge.neighborId[j] = std::numeric_limits<Integer>::max();

                            Integer aFaceId = anElement.faceId(j);
                            CMshFace<Real> aFace = anEdgeMesh.face(aFaceId);

                            if (aFace.hasPair())
                            {
                                Integer aPairFaceId = aFace.pairFaceId();
                                Integer anElementId = anEdgeMesh.face(aPairFaceId).elementId();

                                if (newEdgeIds.find(anElementId) != newEdgeIds.end())
                                {
                                    anAdvEdge.neighborId[j] = newEdgeIds[anElementId];
                                }
                                else
                                    std::cout << "Error: element id = " << anElementId << " not found!" << std::endl;

                            }
                            else
                            {
                                throw std::runtime_error("Boundary is open!");
                            }

                        }

                        m_anAdvFront.push_back(anAdvEdge);

                        this->addEdgeToRtree(anAdvEdge, aTolerance);

                    }

                }

                // Add interior nodes
                for (Integer i = 0; i < static_cast<Integer> (sInteriorPoints.size()); ++i)
                {

                    Integer aNewNodeId = m_surfaceMesh.nextNodeId();
                    CMshNode<Real> aNewNode = sInteriorPoints[i];

                    m_surfaceMesh.addNode(aNewNodeId, aNewNode);

                    SNode anInteriorNode;

                    anInteriorNode.id = i;
                    anInteriorNode.remove = false;

                    anInteriorNode.nodeId = aNewNodeId;

                    m_interiorNodes.push_back(anInteriorNode);

                }

                // Start meshing interior

                Integer maxElem = maxNbElements;

                bool res = true;

                res ? res = this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 0.75, 0.10, false, 0, aTolerance) : res = false;
                res ? res = this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 1.00, 0.05, false, 0, aTolerance) : res = false;
                res ? res = this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 1.50, 0.00, false, 0, aTolerance) : res = false;
                res ? res = this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 2.50, 0.00, false, 0, aTolerance) : res = false;

                m_surfaceMesh.removeDanglingNodes();
                m_surfaceMesh.renumber();

                if (this->frontSize() == 0)
                {
                    m_surfaceMesh.generateFaces(aTolerance);
                    return true;
                }
                else
                {
                    std::cout << "Meshing error!" << std::endl;
                    return false;
                }

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }

        }

        template <typename Real>
        Integer CMshQuadrilateralMesher<Real>::frontSize()
        {

            Integer n = 0;

            for (Integer i = 0; i < static_cast<Integer> (m_anAdvFront.size()); ++i)
            {
                if (!m_anAdvFront[i].remove)
                    n++;
            }

            return n;

        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::advancingFrontQuadMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real sizeFactor, Real shrinkFactor, Real expandFactor, Real minQuality, const bool bCheckDelaunay, Integer firstIndex, const Real aTolerance)
        {

            bool res = false;

            try
            {

                Real x, y;

                meshSizeFunc.removeAllVariables();

                meshSizeFunc.defineVariable("x", x);
                meshSizeFunc.defineVariable("y", y);

                const Real pi = std::acos(-1.0);

                Real sumMeshSize = 0;
                Real averageMeshSize = 0;
                Integer nMeshSize = 0;

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

                        SAdvancingFrontEdge& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                        Integer aNodeId1 = anAdvEdge.nodeId[0];
                        Integer aNodeId2 = anAdvEdge.nodeId[1];

                        CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
                        CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);

                        Integer aNodeId3 = m_anAdvFront[anAdvEdge.neighborId[0]].nodeId[0];
                        Integer aNodeId4 = m_anAdvFront[anAdvEdge.neighborId[1]].nodeId[1];

                        Integer prevPrevEdge = m_anAdvFront[anAdvEdge.neighborId[0]].neighborId[0];
                        Integer nextNextEdge = m_anAdvFront[anAdvEdge.neighborId[1]].neighborId[1];

                        Integer aNodeId5 = m_anAdvFront[prevPrevEdge].nodeId[0];
                        Integer aNodeId6 = m_anAdvFront[nextNextEdge].nodeId[1];

                        if ((aNodeId1 == aNodeId6 || aNodeId2 == aNodeId5) && anAdvEdge.neighborId[0] != nextNextEdge && anAdvEdge.neighborId[1] != prevPrevEdge)
                            continue;

                        CMshNode<Real> aMidNode = (aNode1 + aNode2) * 0.5;

                        x = aMidNode.x();
                        y = aMidNode.y();

                        // Rotate vector by 90º
                        CGeoVector<Real> v = aNode2 - aNode1;

                        Real localMeshSize = std::max(meshSizeFunc.evaluate(), static_cast<Real>(v.norm() * 0.7));

                        sumMeshSize += localMeshSize;
                        nMeshSize++;
                        averageMeshSize = sumMeshSize / nMeshSize;

                        if (averageMeshSize > localMeshSize)
                            localMeshSize = std::min(averageMeshSize, static_cast<Real>(2.0 * localMeshSize));
                        else
                            localMeshSize = std::max(averageMeshSize, static_cast<Real>(0.5 * localMeshSize));

                        v.normalize();

                        CMshNode<Real> aMidNode1 = aMidNode - v * localMeshSize * sizeFactor * 0.5;
                        CMshNode<Real> aMidNode2 = aMidNode + v * localMeshSize * sizeFactor * 0.5;

                        v.rotate(pi*0.5);

                        // Add point to form rectangle with correct spacing
                        CMshNode<Real> aNewNode1 = aMidNode1 + v * localMeshSize * sizeFactor;
                        CMshNode<Real> aNewNode2 = aMidNode2 + v * localMeshSize * sizeFactor;

                        Integer aNewNodeId1 = m_surfaceMesh.nextNodeId() + 0;
                        Integer aNewNodeId2 = m_surfaceMesh.nextNodeId() + 1;

                        // Get closest edges
                        CGeoBoundingBox<Real> aBoundingBox;
                        aBoundingBox.addCoordinate(aNode1);
                        aBoundingBox.addCoordinate(aNode2);
                        aBoundingBox.addCoordinate(aNewNode1);
                        aBoundingBox.addCoordinate(aNewNode2);
                        aBoundingBox.grow(localMeshSize * sizeFactor * 0.5);

                        std::vector<Integer> sEdges;
                        m_tree.find(sEdges, aBoundingBox);

                        // Check if a node exists in proximity
                        std::vector<Integer> sNodes;
                        this->findClosestNodes(sEdges, sNodes);

                        // Meshing priority
                        // Priority = 1: close quadrilateral hole
                        // Priority = 2: close triangular hole
                        // Priority = 3: other nodes in vicinity (3, 4, other node)
                        // Priority = 4: new node forming correct spacing

                        if ((aNodeId3 == aNodeId5 && aNodeId4 == aNodeId6) ||
                            (aNodeId3 == aNodeId6 && aNodeId4 == aNodeId5))
                        {
                            this->addQuadrilateral(anAdvEdge, aNodeId3, aNodeId4, sEdges, aTolerance);
                            res = true;
                            continue;
                        }

                        if (aNodeId3 == aNodeId4)
                        {

                            CMshNode<Real>& aNode3 = m_surfaceMesh.node(aNodeId3);
                            CMshNode<Real>& aNode4 = m_surfaceMesh.node(aNodeId4);

                            CMshTriangle<Real> aNewTriangle;

                            aNewTriangle.addVertex(aNode1);
                            aNewTriangle.addVertex(aNode2);
                            aNewTriangle.addVertex(aNode3);

                            aNewTriangle.calculateArea();

                            if (aNewTriangle.normal().z() > aTolerance)
                            {

                                Integer wNodeId;
                                if (!this->quadrilateralContainsNode(aNode1, aNode2, aNode3, aNode4, wNodeId, sNodes, aTolerance))
                                {

                                    this->addQuadrilateral(anAdvEdge, aNodeId3, aNodeId4, sEdges, aTolerance);
                                    res = true;
                                    continue;

                                }

                            }

                        }

                        Real q1max = 0.0;

                        // If a node exists snap to that node
                        CMshNode<Real> anExistingNode1;
                        Integer anExistingNodeId1 = std::numeric_limits<Integer>::max();

                        for (Integer j = 0; j < static_cast<Integer> (sNodes.size()); ++j)
                        {

                            Integer aNodeId = sNodes[j];

                            if (aNodeId == aNodeId1 ||
                                aNodeId == aNodeId2)
                                continue;

                            CMshNode<Real>& aNode = m_surfaceMesh.node(aNodeId);

                            if (!aBoundingBox.contains(aNode, aTolerance))
                                continue;

                            CMshTriangle<Real> aNewTriangle1;

                            aNewTriangle1.addVertex(aNode1);
                            aNewTriangle1.addVertex(aNode2);
                            aNewTriangle1.addVertex(aNode);

                            aNewTriangle1.calculateArea();

                            if (aNewTriangle1.normal().z() > aTolerance)
                            {

                                CMshQuadrilateral<Real> aNewQuadrilateral1;

                                aNewQuadrilateral1.addVertex(aNode1);
                                aNewQuadrilateral1.addVertex(aNode2);
                                aNewQuadrilateral1.addVertex(aNewNode2);
                                aNewQuadrilateral1.addVertex(aNode);

                                aNewQuadrilateral1.calculateArea();

                                aNewQuadrilateral1.calculateQuality();

                                Real q1 = aNewQuadrilateral1.quality();

                                if (q1 > q1max)
                                    q1 += this->edgeOk(anAdvEdge, aNode1, aNode, sEdges, aTolerance) ? 0.0 : -2.0;

                                Real d = (aNode - aNewNode1).norm();

                                // Use closest node
                                if (d < localMeshSize * sizeFactor * expandFactor * 0.75 && q1 > q1max)
                                {
                                    q1max = q1;
                                    anExistingNodeId1 = aNodeId;
                                    anExistingNode1 = aNode;
                                }

                            }

                        }

                        if (q1max > minQuality)
                        {
                            aNewNode1 = anExistingNode1;
                            aNewNodeId1 = anExistingNodeId1;
                        }

                        Real q2max = 0.0;

                        CMshNode<Real> anExistingNode2;
                        Integer anExistingNodeId2 = std::numeric_limits<Integer>::max();

                        for (Integer j = 0; j < static_cast<Integer> (sNodes.size()); ++j)
                        {

                            Integer aNodeId = sNodes[j];

                            if (aNodeId == aNodeId1 ||
                                aNodeId == aNodeId2)
                                continue;

                            CMshNode<Real>& aNode = m_surfaceMesh.node(aNodeId);

                            if (!aBoundingBox.contains(aNode, aTolerance))
                                continue;

                            CMshTriangle<Real> aNewTriangle2;

                            aNewTriangle2.addVertex(aNode1);
                            aNewTriangle2.addVertex(aNode2);
                            aNewTriangle2.addVertex(aNode);

                            aNewTriangle2.calculateArea();

                            if (aNewTriangle2.normal().z() > aTolerance)
                            {

                                CMshQuadrilateral<Real> aNewQuadrilateral2;

                                aNewQuadrilateral2.addVertex(aNode1);
                                aNewQuadrilateral2.addVertex(aNode2);
                                aNewQuadrilateral2.addVertex(aNode);
                                aNewQuadrilateral2.addVertex(aNewNode1);

                                aNewQuadrilateral2.calculateArea();

                                aNewQuadrilateral2.calculateQuality();

                                Real q2 = aNewQuadrilateral2.quality();

                                if (q2 > q2max)
                                    q2 += this->edgeOk(anAdvEdge, aNode2, aNode, sEdges, aTolerance) ? 0.0 : -2.0;

                                Real d = (aNode - aNewNode2).norm();

                                // Use closest node
                                if (d < localMeshSize * sizeFactor * expandFactor * 0.75 && q2 > q2max)
                                {
                                    q2max = q2;
                                    anExistingNodeId2 = aNodeId;
                                    anExistingNode2 = aNode;
                                }

                            }

                        }

                        if (q2max > minQuality)
                        {
                            aNewNode2 = anExistingNode2;
                            aNewNodeId2 = anExistingNodeId2;
                        }

                        Real q3max = 0.0;

                        CMshQuadrilateral<Real> aNewQuadrilateral3;

                        aNewQuadrilateral3.addVertex(aNode1);
                        aNewQuadrilateral3.addVertex(aNode2);
                        aNewQuadrilateral3.addVertex(aNewNode2);
                        aNewQuadrilateral3.addVertex(aNewNode1);

                        aNewQuadrilateral3.calculateArea();

                        if (aNewQuadrilateral3.normal().z() > aTolerance)
                        {

                            aNewQuadrilateral3.calculateQuality();

                            Real q3 = aNewQuadrilateral3.quality();

                            if (q3 > q3max)
                                q3 += this->edgeOk(anAdvEdge, aNode1, aNewNode1, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q3 > q3max)
                                q3 += this->edgeOk(anAdvEdge, aNode2, aNewNode2, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q3 > q3max)
                                q3 += this->edgeOk(anAdvEdge, aNode1, aNewNode2, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q3 > q3max)
                                q3 += this->edgeOk(anAdvEdge, aNode2, aNewNode1, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q3 > q3max)
                                q3 += this->edgeOk(anAdvEdge, aNewNode1, aNewNode2, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q3 > q3max)
                            {
                                Integer wNodeId;
                                q3 += !this->quadrilateralContainsNode(aNode1, aNode2, aNewNode2, aNewNode1, wNodeId, sNodes, aTolerance) ? 0.0 : -2.0;
                            }

                            q3max = q3;

                        }

                        if (q3max > minQuality)
                        {

                            Integer aNextNodeId = m_surfaceMesh.nextNodeId();

                            if (aNewNodeId1 == aNextNodeId + 0)
                                m_surfaceMesh.addNode(aNewNodeId1, aNewNode1);

                            if (aNewNodeId2 == aNextNodeId + 1)
                                m_surfaceMesh.addNode(aNewNodeId2, aNewNode2);

                            for (Integer k = 0; k < static_cast<Integer> (m_interiorNodes.size()); ++k)
                            {

                                if (m_interiorNodes[k].remove)
                                    continue;

                                if (m_interiorNodes[k].nodeId == aNewNodeId1)
                                    m_interiorNodes[k].remove = true;

                                if (m_interiorNodes[k].nodeId == aNewNodeId2)
                                    m_interiorNodes[k].remove = true;

                            }

                            this->addQuadrilateral(anAdvEdge, aNewNodeId1, aNewNodeId2, sEdges, aTolerance);
                            res = true;

                        }
                        else
                        {

                            CGeoLine<Real> aLine1(aMidNode, aNewNode1);

                            Real dmin1 = findShortestDistance(sEdges, aLine1, anAdvEdgeId, aTolerance);

                            if (dmin1 > localMeshSize * sizeFactor * shrinkFactor * 0.25 && dmin1 < localMeshSize * sizeFactor * expandFactor)
                                aNewNode1 = aNode1 + v * dmin1 * 0.5;

                            CGeoLine<Real> aLine2(aMidNode, aNewNode2);

                            Real dmin2 = findShortestDistance(sEdges, aLine2, anAdvEdgeId, aTolerance);

                            if (dmin2 > localMeshSize * sizeFactor * shrinkFactor * 0.25 && dmin2 < localMeshSize * sizeFactor * expandFactor)
                                aNewNode2 = aNode2 + v * dmin2 * 0.5;

                            CMshQuadrilateral<Real> aNewQuadrilateral4;

                            aNewQuadrilateral4.addVertex(aNode1);
                            aNewQuadrilateral4.addVertex(aNode2);
                            aNewQuadrilateral4.addVertex(aNewNode1);
                            aNewQuadrilateral4.addVertex(aNewNode2);

                            aNewQuadrilateral4.calculateArea();

                            if (aNewQuadrilateral4.normal().z() > aTolerance)
                            {

                                aNewQuadrilateral4.calculateQuality();

                                Real q4 = aNewQuadrilateral4.quality();

                                if (q4 > 1.5 * minQuality)
                                    q4 += this->edgeOk(anAdvEdge, aNode1, aNewNode1, sEdges, aTolerance) ? 0.0 : -2.0;

                                if (q4 > 1.5 * minQuality)
                                    q4 += this->edgeOk(anAdvEdge, aNode2, aNewNode2, sEdges, aTolerance) ? 0.0 : -2.0;

                                if (q4 > 1.5 * minQuality)
                                    q4 += this->edgeOk(anAdvEdge, aNewNode1, aNewNode2, sEdges, aTolerance) ? 0.0 : -2.0;

                                if (q4 > 1.5 * minQuality)
                                {
                                    Integer wNodeId;
                                    q4 += !this->quadrilateralContainsNode(aNode1, aNode2, aNewNode2, aNewNode1, wNodeId, sNodes, aTolerance) ? 0.0 : -2.0;
                                }

                                if (q4 > 1.5 * minQuality)
                                {

                                    // Create a new node
                                    m_surfaceMesh.addNode(aNewNodeId1, aNewNode1);
                                    m_surfaceMesh.addNode(aNewNodeId2, aNewNode2);

                                    this->addQuadrilateral(anAdvEdge, aNewNodeId1, aNewNodeId2, sEdges, aTolerance);
                                    res = true;

                                    if (!m_boundingBox.contains(aNewNode1, aTolerance))
                                    {
                                        throw std::runtime_error("Node is outside boundary!");
                                    }

                                    if (!m_boundingBox.contains(aNewNode2, aTolerance))
                                    {
                                        throw std::runtime_error("Node is outside boundary!");
                                    }

                                }

                            }

                        }

                    }

                }

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
                throw;
            }

            return res;

        }

        template <typename Real>
        CMshMesh<Real>& CMshQuadrilateralMesher<Real>::mesh()
        {

            return m_surfaceMesh;

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::setIntervals(const Integer timeInterval, const Integer dataInterval)
        {

            m_begin = clock();

            m_timeInterval = timeInterval;
            m_dataInterval = dataInterval;

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::stopMeshing()
        {

            m_bStop = true;

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::applyFixedBoundary(CMshMesh<Real>& anEdgeMesh, const Real aTolerance)
        {

            // TODO:

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::flipEdges(const Real aTolerance)
        {

            std::map<Integer, bool> sFlipped;

            for (Integer i = 0; i < m_surfaceMesh.nbElements(); ++i)
            {

                Integer anElementId = m_surfaceMesh.elementId(i);
                sFlipped[anElementId] = false;

            }

            for (Integer i = 0; i < m_surfaceMesh.nbFaces(); ++i)
            {

                Integer aFaceId = m_surfaceMesh.faceId(i);
                CMshFace<Real> aFace = m_surfaceMesh.face(aFaceId);

                if (aFace.faceType() != FT_LINE)
                    continue;

                Integer anElementId = aFace.elementId();

                if (sFlipped[anElementId])
                    continue;

                CMshElement<Real>& anElement = m_surfaceMesh.element(anElementId);

                if (anElement.elementType() != ET_TRIANGLE)
                    continue;

                Integer aNodeId1 = aFace.nodeId(0);
                Integer aNodeId2 = aFace.nodeId(1);

                Integer aNodeId3;

                bool bNodeId3 = false;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {

                    if (anElement.nodeId(j) != aNodeId1 &&
                        anElement.nodeId(j) != aNodeId2)
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
                    CMshFace<Real> aPairFace = m_surfaceMesh.face(pairFaceId);

                    Integer aNeighborId = aPairFace.elementId();

                    if (sFlipped[aNeighborId])
                        continue;

                    CMshElement<Real>& aNeighbor = m_surfaceMesh.element(aNeighborId);

                    if (aNeighbor.elementType() != ET_TRIANGLE)
                        continue;

                    Integer aNodeId4;

                    bool bNodeId4 = false;

                    for (Integer j = 0; j < aNeighbor.nbNodeIds(); ++j)
                    {

                        if (aNeighbor.nodeId(j) != aNodeId1 &&
                            aNeighbor.nodeId(j) != aNodeId2)
                        {
                            aNodeId4 = aNeighbor.nodeId(j);
                            bNodeId4 = true;
                            break;
                        }

                    }

                    if (!bNodeId4)
                        continue;

                    CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
                    CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);
                    CMshNode<Real>& aNode3 = m_surfaceMesh.node(aNodeId3);
                    CMshNode<Real>& aNode4 = m_surfaceMesh.node(aNodeId4);

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
                        m_surfaceMesh.element(anElementId).setNodeId(0, aNodeId3);
                        m_surfaceMesh.element(anElementId).setNodeId(1, aNodeId4);
                        m_surfaceMesh.element(anElementId).setNodeId(2, aNodeId2);

                        m_surfaceMesh.element(aNeighborId).setNodeId(0, aNodeId4);
                        m_surfaceMesh.element(aNeighborId).setNodeId(1, aNodeId3);
                        m_surfaceMesh.element(aNeighborId).setNodeId(2, aNodeId1);

                        sFlipped[anElementId] = true;
                        sFlipped[aNeighborId] = true;

                    }

                }

            }

            m_surfaceMesh.generateFaces(aTolerance);

        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::relaxNodes(const Real aTolerance)
        {

            std::vector<Integer> sElements;

            for (Integer i = 0; i < m_surfaceMesh.nbElements(); ++i)
            {

                Integer anElementId = m_surfaceMesh.elementId(i);
                CMshElement<Real>& anElement = m_surfaceMesh.element(anElementId);

                if (anElement.elementType() != ET_TRIANGLE)
                    continue;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {

                    Integer aPivotNodeId = anElement.nodeId(j);
                    Integer aNextFaceId = anElement.faceId(j);

                    Integer aFirstFaceId = aNextFaceId;

                    CMshNode<Real>& aPivotNode = m_surfaceMesh.node(aPivotNodeId);

                    std::vector<Integer> sNodes;
                    sNodes.push_back(aPivotNodeId);

                    sElements.clear();

                    do
                    {

                        CMshFace<Real>& aNextFace = m_surfaceMesh.face(aNextFaceId);

                        if (!aNextFace.hasPair())
                        {
                            sNodes.clear();
                            sElements.clear();
                            break;
                        }

                        Integer aNextPairFaceId = aNextFace.pairFaceId();
                        CMshFace<Real>& aNextPairFace = m_surfaceMesh.face(aNextPairFaceId);

                        Integer aNextElementId = aNextPairFace.elementId();
                        CMshElement<Real>& aNextElement = m_surfaceMesh.element(aNextElementId);

                        if (aNextElement.elementType() != ET_TRIANGLE)
                            continue;

                        for (Integer k = 0; k < aNextPairFace.nbNodeIds(); ++k)
                        {

                            Integer aNodeId = aNextPairFace.nodeId(k);

                            if (std::find(sNodes.begin(), sNodes.end(), aNodeId) == sNodes.end())
                                sNodes.push_back(aNodeId);

                        }

                        sElements.push_back(aNextElementId);

                        aNextFaceId = aFirstFaceId;

                        for (Integer k = 0; k < aNextElement.nbFaceIds(); ++k)
                        {
                            if (aNextElement.faceId(k) == aNextPairFaceId)
                            {
                                aNextFaceId = aNextElement.faceId((k + 1) % 3);
                                break;
                            }
                        }

                    } while (aNextFaceId != aFirstFaceId);

                    if (sNodes.size() > 1)
                    {

                        CMshNode<Real> aMshNode(0.0, 0.0, 0.0);

                        for (Integer k = 0; k < static_cast<Integer> (sNodes.size()); ++k)
                        {
                            aMshNode += m_surfaceMesh.node(sNodes[k]);
                        }

                        aMshNode /= static_cast<Real>(sNodes.size());

                        Real sumq1 = 0.0;
                        Real sumq2 = 0.0;

                        bool bMoveNode = true;

                        for (Integer k = 0; k < static_cast<Integer> (sElements.size()); ++k)
                        {

                            CMshElement<Real>& aModElement = m_surfaceMesh.element(sElements[k]);

                            Integer aNodeId1 = aModElement.nodeId(0);
                            Integer aNodeId2 = aModElement.nodeId(1);
                            Integer aNodeId3 = aModElement.nodeId(2);

                            CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
                            CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);
                            CMshNode<Real>& aNode3 = m_surfaceMesh.node(aNodeId3);

                            CGeoCoordinate<Real> aVertex1 = aNode1;
                            CGeoCoordinate<Real> aVertex2 = aNode2;
                            CGeoCoordinate<Real> aVertex3 = aNode3;

                            CMshTriangle<Real> aTriangle1;

                            aTriangle1.addVertex(aVertex1);
                            aTriangle1.addVertex(aVertex2);
                            aTriangle1.addVertex(aVertex3);

                            aTriangle1.calculateQuality();

                            Real q1 = aTriangle1.quality();

                            if (aNodeId1 == aPivotNodeId)
                                aVertex1 = aMshNode;

                            if (aNodeId2 == aPivotNodeId)
                                aVertex2 = aMshNode;

                            if (aNodeId3 == aPivotNodeId)
                                aVertex3 = aMshNode;

                            CMshTriangle<Real> aTriangle2;

                            aTriangle2.addVertex(aVertex1);
                            aTriangle2.addVertex(aVertex2);
                            aTriangle2.addVertex(aVertex3);

                            aTriangle2.calculateArea();
                            aTriangle2.calculateQuality();

                            Real q2 = aTriangle2.quality();

                            sumq1 += q1;
                            sumq2 += q2;

                            if (aTriangle2.normal().z() < aTolerance)
                            {
                                bMoveNode = false;
                                break;
                            }

                        }

                        if (bMoveNode && sumq2 > sumq1)
                            aPivotNode = aMshNode;

                    }

                }

            }

        }

    }

}

