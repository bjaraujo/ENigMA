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

#include "GeoCircle.hpp"
#include "GeoLine.hpp"
#include "MshQuadrilateral.hpp"
#include "MshTriangle.hpp"

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        CMshQuadrilateralMesher<Real>::CMshQuadrilateralMesher()
            : CMshTriangleMesher<Real>()
        {
        }

        template <typename Real>
        CMshQuadrilateralMesher<Real>::~CMshQuadrilateralMesher()
        {
        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::cleanDuplicateEdges(std::vector<Integer>& sEdges, const Real aTolerance)
        {
            if (this->m_nextEdgeId < 3)
                return;

            bool bAdjustConnectivity = false;

            // Check last three edges
            for (Integer i = 0; i < 3; ++i)
            {
                Integer anAdvEdgeId = this->m_nextEdgeId - i - 1;

                SMshAdvancingFrontEdge<Real>& anAdvEdge = this->m_anAdvFront[anAdvEdgeId];

                if (anAdvEdge.remove)
                    continue;

                Integer aDuplicateEdgeId;

                if (this->edgeExists(anAdvEdge, aDuplicateEdgeId, sEdges))
                {
                    SMshAdvancingFrontEdge<Real>& aDuplicateEdge = this->m_anAdvFront[aDuplicateEdgeId];

                    this->removeEdge(anAdvEdge, aTolerance);
                    this->removeEdge(aDuplicateEdge, aTolerance);

                    bAdjustConnectivity = true;
                }
            }

            if (bAdjustConnectivity)
                this->adjustConnectivity(sEdges);
        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::quadrilateralContainsNode(CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, CMshNode<Real>& aNode4, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance)
        {
            CMshQuadrilateral<Real> aQuadrilateral;
            aQuadrilateral.reset();
            aQuadrilateral.addVertex(aNode1);
            aQuadrilateral.addVertex(aNode2);
            aQuadrilateral.addVertex(aNode3);
            aQuadrilateral.addVertex(aNode4);

            for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j)
            {
                aNodeId = sNodes[j];

                CMshNode<Real>& aNode = this->m_surfaceMesh.node(aNodeId);

                if ((aNode - aNode1).norm() < aTolerance || (aNode - aNode2).norm() < aTolerance || (aNode - aNode3).norm() < aTolerance || (aNode - aNode4).norm() < aTolerance)
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
        void CMshQuadrilateralMesher<Real>::addQuadrilateral(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Integer aNodeId3, const Integer aNodeId4, std::vector<Integer>& sEdges, const Real aTolerance)
        {
            Integer aNodeId1 = anAdvEdge.nodeId[0];
            Integer aNodeId2 = anAdvEdge.nodeId[1];

            if (aNodeId3 == aNodeId1 || aNodeId3 == aNodeId2 || aNodeId4 == aNodeId1 || aNodeId4 == aNodeId2)
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

                Integer aNewElementId = this->m_surfaceMesh.nextElementId();

                this->m_surfaceMesh.addElement(aNewElementId, aNewElement);

                SMshAdvancingFrontEdge<Real>& aPrevEdge = this->m_anAdvFront[anAdvEdge.neighborId[0]];
                SMshAdvancingFrontEdge<Real>& aNextEdge = this->m_anAdvFront[anAdvEdge.neighborId[1]];

                Integer aNewEdgeId1 = this->m_nextEdgeId++;
                Integer aNewEdgeId2 = this->m_nextEdgeId++;

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

                aNewEdge1.build(this->m_surfaceMesh);

                // Correct connectivity
                aPrevEdge.neighborId[1] = aNewEdgeId1;

                // Add edge to rtree
                this->addEdgeToRtree(aNewEdge1, aTolerance);

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

                aNewEdge2.build(this->m_surfaceMesh);

                // Correct connectivity
                aNextEdge.neighborId[0] = aNewEdgeId2;

                // Add edge to rtree
                this->addEdgeToRtree(aNewEdge2, aTolerance);

                // Remove this edge
                this->removeEdge(anAdvEdge, aTolerance);

                this->m_anAdvFront.push_back(aNewEdge1);
                this->m_anAdvFront.push_back(aNewEdge2);

                this->cleanDuplicateEdges(sEdges, aTolerance);
            }
            else
            {
                // Check if convex
                std::vector<Integer> sNodes = { aNodeId1, aNodeId2, aNodeId4, aNodeId3 };

                Integer wIndex = -1;
                for (Integer i = 0; i < 4; i++)
                {

                    CGeoVector<Real> v1 = this->m_surfaceMesh.node(sNodes[(i + 1) % 4]) - this->m_surfaceMesh.node(sNodes[(i + 0) % 4]);
                    CGeoVector<Real> v2 = this->m_surfaceMesh.node(sNodes[(i + 2) % 4]) - this->m_surfaceMesh.node(sNodes[(i + 1) % 4]);

                    Real z = v1.cross(v2).z();

                    if (z <= 0.0)
                    {
                        wIndex = i;
                        break;
                    }
                }

                Integer aNewElementId = -1;

                if (wIndex == -1)
                {
                    // Add new quadrilateral
                    CMshElement<Real> aNewElement(ET_QUADRILATERAL);
                    aNewElement.addNodeId(aNodeId1);
                    aNewElement.addNodeId(aNodeId2);
                    aNewElement.addNodeId(aNodeId4);
                    aNewElement.addNodeId(aNodeId3);

                    aNewElementId = this->m_surfaceMesh.nextElementId();

                    this->m_surfaceMesh.addElement(aNewElementId, aNewElement);
                }
                else
                {
                    // Add new triangle
                    CMshElement<Real> aNewElement1(ET_TRIANGLE);
                    aNewElement1.addNodeId(sNodes[(wIndex + 0) % 4]); // aNodeId1
                    aNewElement1.addNodeId(sNodes[(wIndex + 1) % 4]); // aNodeId2
                    aNewElement1.addNodeId(sNodes[(wIndex + 3) % 4]); // aNodeId3

                    aNewElementId = this->m_surfaceMesh.nextElementId();

                    this->m_surfaceMesh.addElement(aNewElementId, aNewElement1);

                    // Add new triangle
                    CMshElement<Real> aNewElement2(ET_TRIANGLE);
                    aNewElement2.addNodeId(sNodes[(wIndex + 1) % 4]); // aNodeId2
                    aNewElement2.addNodeId(sNodes[(wIndex + 2) % 4]); // aNodeId4
                    aNewElement2.addNodeId(sNodes[(wIndex + 3) % 4]); // aNodeId3

                    aNewElementId = this->m_surfaceMesh.nextElementId();

                    this->m_surfaceMesh.addElement(aNewElementId, aNewElement2);
                }

                SMshAdvancingFrontEdge<Real>& aPrevEdge = this->m_anAdvFront[anAdvEdge.neighborId[0]];
                SMshAdvancingFrontEdge<Real>& aNextEdge = this->m_anAdvFront[anAdvEdge.neighborId[1]];

                Integer aNewEdgeId1 = this->m_nextEdgeId++;
                Integer aNewEdgeId2 = this->m_nextEdgeId++;
                Integer aNewEdgeId3 = this->m_nextEdgeId++;

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

                aNewEdge1.build(this->m_surfaceMesh);

                // Correct connectivity
                aPrevEdge.neighborId[1] = aNewEdgeId1;

                // Add edge to rtree
                this->addEdgeToRtree(aNewEdge1, aTolerance);

                // Add edge 2
                SMshAdvancingFrontEdge<Real> aNewEdge2;
                aNewEdge2.id = aNewEdgeId2;
                aNewEdge2.remove = false;
                aNewEdge2.boundary = false;
                aNewEdge2.nodeId[0] = aNodeId3;
                aNewEdge2.nodeId[1] = aNodeId4;
                aNewEdge2.neighborId[0] = aNewEdgeId1;
                aNewEdge2.neighborId[1] = aNewEdgeId3;
                aNewEdge2.elementId = aNewElementId;

                aNewEdge2.build(this->m_surfaceMesh);

                // Add edge to rtree
                this->addEdgeToRtree(aNewEdge2, aTolerance);

                // Add edge 3
                SMshAdvancingFrontEdge<Real> aNewEdge3;
                aNewEdge3.id = aNewEdgeId3;
                aNewEdge3.remove = false;
                aNewEdge3.boundary = false;
                aNewEdge3.nodeId[0] = aNodeId4;
                aNewEdge3.nodeId[1] = aNodeId2;
                aNewEdge3.neighborId[0] = aNewEdgeId2;
                aNewEdge3.neighborId[1] = anAdvEdge.neighborId[1];
                aNewEdge3.elementId = aNewElementId;

                aNewEdge3.build(this->m_surfaceMesh);

                // Correct connectivity
                aNextEdge.neighborId[0] = aNewEdgeId3;

                // Add edge to rtree
                this->addEdgeToRtree(aNewEdge3, aTolerance);

                // Remove this edge
                this->removeEdge(anAdvEdge, aTolerance);

                this->m_anAdvFront.push_back(aNewEdge1);
                this->m_anAdvFront.push_back(aNewEdge2);
                this->m_anAdvFront.push_back(aNewEdge3);

                this->cleanDuplicateEdges(sEdges, aTolerance);
            }
        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::generate(const CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, Real meshSize, Real minMeshSize, Real maxMeshSize, const Real aTolerance)
        {
            std::vector<CGeoCoordinate<Real>> sInteriorPoints;

            ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
            aAnaFunction.set(meshSize);

            return this->generate(anEdgeMesh, maxNbElements, sInteriorPoints, aAnaFunction, minMeshSize, maxMeshSize, aTolerance);
        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::generate(const CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, Real meshSize, Real minMeshSize, Real maxMeshSize, const Real aTolerance)
        {
            ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
            aAnaFunction.set(meshSize);

            return this->generate(anEdgeMesh, maxNbElements, sInteriorPoints, aAnaFunction, minMeshSize, maxMeshSize, aTolerance);
        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::generate(const CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minMeshSize, Real maxMeshSize, const Real aTolerance)
        {
            this->m_previousNbElements = 0;

            this->m_bStop = false;

            // Add boundary nodes to surface mesh
            this->m_surfaceMesh.reset();
            this->m_boundingBox.reset();

            for (Integer i = 0; i < anEdgeMesh.nbNodes(); ++i)
            {
                Integer aNodeId = anEdgeMesh.nodeId(i);
                CMshNode<Real> aNode = anEdgeMesh.node(aNodeId);

                this->m_surfaceMesh.addNode(aNodeId, aNode);

                // Add to bounding box
                this->m_boundingBox.addCoordinate(aNode);
            }

            this->m_boundingBox.grow(aTolerance);

            // Add boundary to advancing front and rtree
            this->m_anAdvFront.clear();

            this->m_anAdvFront.reserve(anEdgeMesh.nbElements() * 50);

            this->m_tree.reset();

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

            this->m_nextEdgeId = 0;

            for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i)
            {
                Integer anAdvEdgeId = anEdgeMesh.elementId(i);

                const CMshElement<Real>& anElement = anEdgeMesh.element(anAdvEdgeId);

                if (anElement.elementType() == ET_BEAM)
                {
                    SMshAdvancingFrontEdge<Real> anAdvEdge;

                    anAdvEdge.id = this->m_nextEdgeId++;
                    anAdvEdge.remove = false;
                    anAdvEdge.boundary = true;

                    for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                        anAdvEdge.nodeId[j] = anElement.nodeId(j);

                    anAdvEdge.elementId = std::numeric_limits<Integer>::max();

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

                    this->m_anAdvFront.push_back(anAdvEdge);

                    this->addEdgeToRtree(anAdvEdge, aTolerance);
                }
            }

            // Add interior nodes
            for (Integer i = 0; i < static_cast<Integer>(sInteriorPoints.size()); ++i)
            {
                Integer aNewNodeId = this->m_surfaceMesh.nextNodeId();
                CMshNode<Real> aNewNode = sInteriorPoints[i];

                this->m_surfaceMesh.addNode(aNewNodeId, aNewNode);

                SNode<Real> anInteriorNode;

                anInteriorNode.id = i;
                anInteriorNode.remove = false;

                anInteriorNode.nodeId = aNewNodeId;

                this->m_interiorNodes.push_back(anInteriorNode);
            }

            // Start meshing interior

            Integer maxElem = maxNbElements;

            bool res = true;

            this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 0.75, 0.05, false, false, aTolerance);
            this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 0.75, 0.05, true, false, aTolerance);
            this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 1.00, 0.02, true, false, aTolerance);
            this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 1.50, 0.00, true, false, aTolerance);
            this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 2.50, 0.00, true, false, aTolerance);
            this->advancingFrontQuadMeshing(meshSizeFunc, maxElem, minMeshSize, maxMeshSize, 1.00, 1.00, 3.50, 0.00, true, false, aTolerance);

            if (this->frontSize() > 0)
            {
                std::cout << "Meshing error!" << std::endl;
                throw(this->m_anAdvFront);
            }

            this->m_surfaceMesh.removeDanglingNodes();
            this->m_surfaceMesh.renumber();
            this->m_surfaceMesh.generateFaces(aTolerance);

            return true;
        }

        template <typename Real>
        bool CMshQuadrilateralMesher<Real>::advancingFrontQuadMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real minMeshSize, Real maxMeshSize, Real sizeFactor, Real shrinkFactor, Real expandFactor, Real minQuality, const bool bAddNodes, const bool bCheckDelaunay, const Real aTolerance)
        {
            bool res = false;

            Real x, y;

            meshSizeFunc.removeAllVariables();

            meshSizeFunc.defineVariable("x", x);
            meshSizeFunc.defineVariable("y", y);

            static const Real pi = std::acos(-1.0);

            std::vector<Integer> sNodes;
            std::vector<Integer> sEdges;

            CGeoVector<Real> a;
            CGeoVector<Real> v;

            CGeoLine<Real> aLine1;

            CMshNode<Real> aMidNode;
            CMshNode<Real> aMidNode1;
            CMshNode<Real> aMidNode2;
            CMshNode<Real> aNewNode1;
            CMshNode<Real> aNewNode2;
            CMshNode<Real> anExistingNode1;
            CMshNode<Real> anExistingNode2;

            CGeoBoundingBox<Real> aBoundingBox;

            CMshTriangle<Real> aNewTriangle;
            CMshTriangle<Real> aNewTriangle1;
            CMshTriangle<Real> aNewTriangle2;

            CMshQuadrilateral<Real> aNewQuadrilateral1;
            CMshQuadrilateral<Real> aNewQuadrilateral2;
            CMshQuadrilateral<Real> aNewQuadrilateral3;
            CMshQuadrilateral<Real> aNewQuadrilateral4;

            for (Integer i = 0; i < static_cast<Integer>(this->m_anAdvFront.size()); ++i)
            {
                this->checkUpdate();

                if (this->m_bStop)
                    return false;

                if (maxNbElements > 0)
                {
                    if (this->m_surfaceMesh.nbElements() >= maxNbElements)
                    {
                        std::cout << "Max number of elements (" << maxNbElements << ") reached!" << std::endl;
                        return false;
                    }
                }

                if (!this->m_anAdvFront[i].remove)
                {
                    Integer anAdvEdgeId = i;

                    SMshAdvancingFrontEdge<Real>& anAdvEdge = this->m_anAdvFront[anAdvEdgeId];

                    Integer aNodeId1 = anAdvEdge.nodeId[0];
                    Integer aNodeId2 = anAdvEdge.nodeId[1];

                    CMshNode<Real>& aNode1 = this->m_surfaceMesh.node(aNodeId1);
                    CMshNode<Real>& aNode2 = this->m_surfaceMesh.node(aNodeId2);

                    Integer aNodeId3 = this->m_anAdvFront[anAdvEdge.neighborId[0]].nodeId[0];
                    Integer aNodeId4 = this->m_anAdvFront[anAdvEdge.neighborId[1]].nodeId[1];

                    Integer prevPrevEdge = this->m_anAdvFront[anAdvEdge.neighborId[0]].neighborId[0];
                    Integer nextNextEdge = this->m_anAdvFront[anAdvEdge.neighborId[1]].neighborId[1];

                    Integer aNodeId5 = this->m_anAdvFront[prevPrevEdge].nodeId[0];
                    Integer aNodeId6 = this->m_anAdvFront[nextNextEdge].nodeId[1];

                    if ((aNodeId1 == aNodeId6 || aNodeId2 == aNodeId5) && anAdvEdge.neighborId[0] != nextNextEdge && anAdvEdge.neighborId[1] != prevPrevEdge)
                        continue;

                    aMidNode = (aNode1 + aNode2) * 0.5;

                    x = aMidNode.x();
                    y = aMidNode.y();

                    a = aNode2 - aNode1;

                    Real requiredMeshSize = meshSizeFunc.evaluate();
                    Real localMeshSize = static_cast<Real>(a.norm());

                    Real aFactor = std::max<Real>(std::min<Real>(requiredMeshSize / (localMeshSize + aTolerance * aTolerance), 2.0), 0.5);

                    localMeshSize *= aFactor;
                    localMeshSize = std::max(localMeshSize, minMeshSize);
                    localMeshSize = std::min(localMeshSize, maxMeshSize);

                    Real baseHeightSize = localMeshSize; // Equilateral rectangle (height to edge ratio)

                    v = a;
                    v.normalize();

                    aMidNode1 = aMidNode - v * baseHeightSize * sizeFactor * 0.5;
                    aMidNode2 = aMidNode + v * baseHeightSize * sizeFactor * 0.5;

                    // Rotate vector by 90 degrees
                    v.rotate(pi * 0.5);

                    // Add point to form rectangle with correct spacing
                    aNewNode1 = aMidNode1 + v * baseHeightSize * sizeFactor;
                    aNewNode2 = aMidNode2 + v * baseHeightSize * sizeFactor;

                    Integer aNewNodeId1 = this->m_surfaceMesh.nextNodeId() + 0;
                    Integer aNewNodeId2 = this->m_surfaceMesh.nextNodeId() + 1;

                    // Get closest edges
                    aBoundingBox.reset();
                    aBoundingBox.addCoordinate(aNode1);
                    aBoundingBox.addCoordinate(aNode2);
                    aBoundingBox.addCoordinate(aNewNode1);
                    aBoundingBox.addCoordinate(aNewNode2);
                    aBoundingBox.grow(baseHeightSize * sizeFactor * 0.5);

                    sEdges.clear();
                    this->m_tree.find(sEdges, aBoundingBox);

                    // Check if a node exists in proximity
                    this->findClosestNodes(sEdges, sNodes);

                    // Meshing priority
                    // Priority = 1: close quadrilateral hole
                    // Priority = 2: close triangular hole
                    // Priority = 3: other nodes in vicinity (3, 4, other node)
                    // Priority = 4: new node forming correct spacing

                    if ((aNodeId3 == aNodeId5 && aNodeId4 == aNodeId6) || (aNodeId3 == aNodeId6 && aNodeId4 == aNodeId5))
                    {
                        this->addQuadrilateral(anAdvEdge, aNodeId3, aNodeId4, sEdges, aTolerance);
                        res = true;
                        continue;
                    }

                    if (aNodeId3 == aNodeId4)
                    {
                        CMshNode<Real>& aNode3 = this->m_surfaceMesh.node(aNodeId3);
                        CMshNode<Real>& aNode4 = this->m_surfaceMesh.node(aNodeId4);

                        aNewTriangle.reset();
                        aNewTriangle.addVertex(aNode1);
                        aNewTriangle.addVertex(aNode2);
                        aNewTriangle.addVertex(aNode3);

                        aNewTriangle.calculateArea();

                        if (aNewTriangle.normal().z() > aTolerance)
                        {
                            Integer wNodeId;

                            Real q0 = 1.0;

                            if (q0 > 0.0)
                                q0 += this->edgeOk(anAdvEdge, aNode1, aNode3, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q0 > 0.0)
                                q0 += this->edgeOk(anAdvEdge, aNode2, aNode4, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (!this->quadrilateralContainsNode(aNode1, aNode2, aNode3, aNode4, wNodeId, sNodes, aTolerance) && q0 > 0.0)
                            {
                                this->addQuadrilateral(anAdvEdge, aNodeId3, aNodeId4, sEdges, aTolerance);
                                res = true;
                                continue;
                            }
                        }
                    }

                    Real q1max = 0.0;

                    // If a node exists snap to that node
                    Integer anExistingNodeId1 = std::numeric_limits<Integer>::max();

                    for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j)
                    {
                        Integer aNodeId = sNodes[j];

                        if (aNodeId == aNodeId1 || aNodeId == aNodeId2)
                            continue;

                        CMshNode<Real>& aNode = this->m_surfaceMesh.node(aNodeId);

                        if (!aBoundingBox.contains(aNode, aTolerance))
                            continue;

                        aNewTriangle1.reset();
                        aNewTriangle1.addVertex(aNode1);
                        aNewTriangle1.addVertex(aNode2);
                        aNewTriangle1.addVertex(aNode);

                        aNewTriangle1.calculateArea();

                        if (aNewTriangle1.normal().z() > aTolerance)
                        {
                            aNewQuadrilateral1.reset();
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
                            if (d < baseHeightSize * sizeFactor * expandFactor * 0.75 && q1 > q1max)
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

                    Integer anExistingNodeId2 = std::numeric_limits<Integer>::max();

                    for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j)
                    {
                        Integer aNodeId = sNodes[j];

                        if (aNodeId == aNodeId1 || aNodeId == aNodeId2)
                            continue;

                        CMshNode<Real>& aNode = this->m_surfaceMesh.node(aNodeId);

                        if (!aBoundingBox.contains(aNode, aTolerance))
                            continue;

                        aNewTriangle2.reset();
                        aNewTriangle2.addVertex(aNode1);
                        aNewTriangle2.addVertex(aNode2);
                        aNewTriangle2.addVertex(aNode);

                        aNewTriangle2.calculateArea();

                        if (aNewTriangle2.normal().z() > aTolerance)
                        {
                            aNewQuadrilateral2.reset();
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
                            if (d < baseHeightSize * sizeFactor * expandFactor * 0.75 && q2 > q2max)
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

                    aNewQuadrilateral3.reset();
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
                        Integer aNextNodeId = this->m_surfaceMesh.nextNodeId();

                        if (aNewNodeId1 == aNextNodeId + 0)
                            this->m_surfaceMesh.addNode(aNewNodeId1, aNewNode1);

                        if (aNewNodeId2 == aNextNodeId + 1)
                            this->m_surfaceMesh.addNode(aNewNodeId2, aNewNode2);

                        for (Integer k = 0; k < static_cast<Integer>(this->m_interiorNodes.size()); ++k)
                        {
                            if (this->m_interiorNodes[k].remove)
                                continue;

                            if (this->m_interiorNodes[k].nodeId == aNewNodeId1)
                                this->m_interiorNodes[k].remove = true;

                            if (this->m_interiorNodes[k].nodeId == aNewNodeId2)
                                this->m_interiorNodes[k].remove = true;
                        }

                        this->addQuadrilateral(anAdvEdge, aNewNodeId1, aNewNodeId2, sEdges, aTolerance);
                        res = true;
                    }
                    else if (bAddNodes)
                    {
                        aLine1.reset();
                        aLine1.setStartPoint(aMidNode);
                        aLine1.setEndPoint(aNewNode1);

                        Real dmin1 = this->findShortestDistance(sEdges, aLine1, anAdvEdgeId, aTolerance);

                        if (dmin1 > baseHeightSize * sizeFactor * shrinkFactor * 0.25 && dmin1 < baseHeightSize * sizeFactor * expandFactor)
                            aNewNode1 = aNode1 + v * dmin1 * 0.5;

                        CGeoLine<Real> aLine2(aMidNode, aNewNode2);

                        Real dmin2 = this->findShortestDistance(sEdges, aLine2, anAdvEdgeId, aTolerance);

                        if (dmin2 > baseHeightSize * sizeFactor * shrinkFactor * 0.25 && dmin2 < baseHeightSize * sizeFactor * expandFactor)
                            aNewNode2 = aNode2 + v * dmin2 * 0.5;

                        aNewQuadrilateral4.reset();
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
                                this->m_surfaceMesh.addNode(aNewNodeId1, aNewNode1);
                                this->m_surfaceMesh.addNode(aNewNodeId2, aNewNode2);

                                this->addQuadrilateral(anAdvEdge, aNewNodeId1, aNewNodeId2, sEdges, aTolerance);
                                res = true;

                                if (!this->m_boundingBox.contains(aNewNode1, aTolerance))
                                {
                                    throw std::runtime_error("Node is outside boundary!");
                                }

                                if (!this->m_boundingBox.contains(aNewNode2, aTolerance))
                                {
                                    throw std::runtime_error("Node is outside boundary!");
                                }
                            }
                        }
                    }
                }
            }

            return res;
        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::flipEdges(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aTolerance)
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

                if (sFlipped[anElementId])
                    continue;

                CMshElement<Real>& anElement = aMesh.element(anElementId);

                if (anElement.elementType() != ET_TRIANGLE)
                    continue;

                Integer aNodeId1 = aFace.nodeId(0);
                Integer aNodeId2 = aFace.nodeId(1);
                Integer aNodeId3;
                Integer aNodeId4;
                Integer aNodeId5;

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
                    const CMshFace<Real> aPairFace = aMesh.face(pairFaceId);

                    Integer aNeighborId = aPairFace.elementId();

                    if (sFlipped[aNeighborId])
                        continue;

                    CMshElement<Real>& aNeighbor = aMesh.element(aNeighborId);

                    if (aNeighbor.elementType() != ET_QUADRILATERAL)
                        continue;

                    Integer wIndex = -1;
                    for (Integer j = 0; j < aNeighbor.nbNodeIds(); ++j)
                    {
                        if (aNeighbor.nodeId((j + 1) % 4) == aNodeId1 && aNeighbor.nodeId((j + 0) % 4) == aNodeId2)
                        {
                            wIndex = j;
                            break;
                        }
                    }

                    if (wIndex == -1)
                    {
                        continue;
                    }

                    aNodeId4 = aNeighbor.nodeId((wIndex + 3) % 4);
                    aNodeId5 = aNeighbor.nodeId((wIndex + 2) % 4);

                    const CMshNode<Real>& aNode1 = aMesh.node(aNodeId1);
                    const CMshNode<Real>& aNode2 = aMesh.node(aNodeId2);
                    const CMshNode<Real>& aNode3 = aMesh.node(aNodeId3);
                    const CMshNode<Real>& aNode4 = aMesh.node(aNodeId4);
                    const CMshNode<Real>& aNode5 = aMesh.node(aNodeId5);

                    // Original
                    CMshTriangle<Real> aTriangle1;
                    aTriangle1.addVertex(aNode1);
                    aTriangle1.addVertex(aNode2);
                    aTriangle1.addVertex(aNode3);

                    aTriangle1.calculateArea();
                    aTriangle1.calculateQuality();
                    Real q1 = aTriangle1.quality();

                    CMshQuadrilateral<Real> aQuadrilateral1;
                    aQuadrilateral1.addVertex(aNode2);
                    aQuadrilateral1.addVertex(aNode1);
                    aQuadrilateral1.addVertex(aNode5);
                    aQuadrilateral1.addVertex(aNode4);

                    aQuadrilateral1.calculateQuality();
                    Real q2 = aQuadrilateral1.quality();

                    // Modified 1
                    CMshTriangle<Real> aTriangle2;
                    aTriangle2.addVertex(aNode2);
                    aTriangle2.addVertex(aNode3);
                    aTriangle2.addVertex(aNode4);

                    aTriangle2.calculateArea();
                    aTriangle2.calculateQuality();
                    Real q3 = aTriangle2.quality();

                    CMshQuadrilateral<Real> aQuadrilateral2;
                    aQuadrilateral2.addVertex(aNode1);
                    aQuadrilateral2.addVertex(aNode5);
                    aQuadrilateral2.addVertex(aNode4);
                    aQuadrilateral2.addVertex(aNode3);

                    aQuadrilateral2.calculateQuality();
                    Real q4 = aQuadrilateral2.quality();

                    if (std::min(q3, q4) > std::min(q1, q2) && aTriangle2.normal().z() > aTolerance && aQuadrilateral2.normal().z() > aTolerance)
                    {

                        // Check if convex
                        std::vector<Integer> sNodes = { aNodeId1, aNodeId5, aNodeId4, aNodeId3 };

                        bool isConvex = true;
                        for (Integer i = 0; i < 4; i++)
                        {

                            CGeoVector<Real> v1 = aMesh.node(sNodes[(i + 1) % 4]) - aMesh.node(sNodes[(i + 0) % 4]);
                            CGeoVector<Real> v2 = aMesh.node(sNodes[(i + 2) % 4]) - aMesh.node(sNodes[(i + 1) % 4]);

                            Real z = v1.cross(v2).z();

                            if (z <= 0.0)
                            {
                                isConvex = false;
                                break;
                            }
                        }

                        if (isConvex)
                        {
                            // Do flip
                            aMesh.element(anElementId).setNodeId(0, aNodeId2);
                            aMesh.element(anElementId).setNodeId(1, aNodeId3);
                            aMesh.element(anElementId).setNodeId(2, aNodeId4);

                            aMesh.element(aNeighborId).setNodeId(0, aNodeId1);
                            aMesh.element(aNeighborId).setNodeId(1, aNodeId5);
                            aMesh.element(aNeighborId).setNodeId(2, aNodeId4);
                            aMesh.element(aNeighborId).setNodeId(3, aNodeId3);

                            sFlipped[anElementId] = true;
                            sFlipped[aNeighborId] = true;
                            continue;
                        }
                    }

                    // Modified 2
                    CMshTriangle<Real> aTriangle3;
                    aTriangle3.addVertex(aNode1);
                    aTriangle3.addVertex(aNode5);
                    aTriangle3.addVertex(aNode3);

                    aTriangle3.calculateArea();
                    aTriangle3.calculateQuality();
                    Real q5 = aTriangle3.quality();

                    CMshQuadrilateral<Real> aQuadrilateral3;
                    aQuadrilateral3.addVertex(aNode3);
                    aQuadrilateral3.addVertex(aNode5);
                    aQuadrilateral3.addVertex(aNode4);
                    aQuadrilateral3.addVertex(aNode2);

                    aQuadrilateral3.calculateQuality();
                    Real q6 = aQuadrilateral3.quality();

                    if (std::min(q6, q5) > std::min(q1, q2) && aTriangle3.normal().z() > aTolerance && aQuadrilateral3.normal().z() > aTolerance)
                    {

                        // Check if convex
                        std::vector<Integer> sNodes = { aNodeId3, aNodeId5, aNodeId4, aNodeId2 };

                        bool isConvex = true;
                        for (Integer i = 0; i < 4; i++)
                        {

                            CGeoVector<Real> v1 = aMesh.node(sNodes[(i + 1) % 4]) - aMesh.node(sNodes[(i + 0) % 4]);
                            CGeoVector<Real> v2 = aMesh.node(sNodes[(i + 2) % 4]) - aMesh.node(sNodes[(i + 1) % 4]);

                            Real z = v1.cross(v2).z();

                            if (z <= 0.0)
                            {
                                isConvex = false;
                                break;
                            }
                        }

                        if (isConvex)
                        {
                            // Do flip
                            aMesh.element(anElementId).setNodeId(0, aNodeId1);
                            aMesh.element(anElementId).setNodeId(1, aNodeId5);
                            aMesh.element(anElementId).setNodeId(2, aNodeId3);

                            aMesh.element(aNeighborId).setNodeId(0, aNodeId3);
                            aMesh.element(aNeighborId).setNodeId(1, aNodeId5);
                            aMesh.element(aNeighborId).setNodeId(2, aNodeId4);
                            aMesh.element(aNeighborId).setNodeId(3, aNodeId2);

                            sFlipped[anElementId] = true;
                            sFlipped[aNeighborId] = true;
                            continue;
                        }
                    }
                }
            }

            aMesh.generateFaces(aTolerance);
        }

        template <typename Real>
        void CMshQuadrilateralMesher<Real>::relaxNodes(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aTolerance)
        {
            std::vector<Integer> sElements;

            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {
                Integer anElementId = aMesh.elementId(i);
                const CMshElement<Real>& anElement = aMesh.element(anElementId);

                if (anElement.elementType() != ET_TRIANGLE && anElement.elementType() != ET_QUADRILATERAL)
                    continue;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    Integer aPivotNodeId = anElement.nodeId(j);
                    Integer aNextFaceId = anElement.faceId(j);

                    Integer aFirstFaceId = aNextFaceId;

                    CMshNode<Real>& aPivotNode = aMesh.node(aPivotNodeId);

                    std::vector<Integer> sNodes;
                    sNodes.push_back(aPivotNodeId);

                    sElements.clear();

                    Integer iter = 0;

                    do
                    {
                        const CMshFace<Real>& aNextFace = aMesh.face(aNextFaceId);

                        if (!aNextFace.hasPair())
                        {
                            sNodes.clear();
                            sElements.clear();
                            break;
                        }

                        Integer aNextPairFaceId = aNextFace.pairFaceId();
                        const CMshFace<Real>& aNextPairFace = aMesh.face(aNextPairFaceId);

                        Integer aNextElementId = aNextPairFace.elementId();
                        const CMshElement<Real>& aNextElement = aMesh.element(aNextElementId);

                        if (aNextElement.elementType() != ET_TRIANGLE && aNextElement.elementType() != ET_QUADRILATERAL)
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

                        iter++;

                        if (iter > 1000)
                            return;

                    } while (aNextFaceId != aFirstFaceId);

                    if (sNodes.size() > 1)
                    {
                        CMshNode<Real> aMshNode(0.0, 0.0, 0.0);

                        for (Integer k = 0; k < static_cast<Integer>(sNodes.size()); ++k)
                        {
                            aMshNode += aMesh.node(sNodes[k]);
                        }

                        aMshNode /= static_cast<Real>(sNodes.size());

                        Real sumq1 = 0.0;
                        Real sumq2 = 0.0;

                        bool bMoveNode = true;

                        for (Integer k = 0; k < static_cast<Integer>(sElements.size()); ++k)
                        {
                            CMshElement<Real>& aModElement = aMesh.element(sElements[k]);

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
