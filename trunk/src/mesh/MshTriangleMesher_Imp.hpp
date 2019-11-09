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
#include "MshMeshQuery.hpp"
#include "MshTriangle.hpp"

using namespace ENigMA::geometry;

namespace ENigMA {

namespace mesh {

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

        if (onUpdate != nullptr) {
            m_end = clock();

            if (static_cast<Real>(m_end - m_begin) / CLOCKS_PER_SEC > m_timeInterval) {
                m_begin = m_end;

                if (m_surfaceMesh.nbElements() > m_previousNbElements + m_dataInterval) {
                    m_previousNbElements = m_surfaceMesh.nbElements();
                    onUpdate(true);
                } else
                    onUpdate(false);
            }
        }
    }

    template <typename Real>
    void CMshTriangleMesher<Real>::removeEdge(SAdvancingFrontEdge& anAdvEdge, const Real aTolerance)
    {

        anAdvEdge.remove = true;
        removeEdgeFromRtree(anAdvEdge, aTolerance);
    }

    template <typename Real>
    void CMshTriangleMesher<Real>::addEdgeToRtree(SAdvancingFrontEdge& anAdvEdge, const Real aTolerance)
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
    void CMshTriangleMesher<Real>::removeEdgeFromRtree(SAdvancingFrontEdge& anAdvEdge, const Real aTolerance)
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
    bool CMshTriangleMesher<Real>::edgeExists(SAdvancingFrontEdge& anAdvEdge, Integer& aDuplicateEdgeId, std::vector<Integer>& sEdges)
    {

        std::vector<Integer> sNodeIds;

        sNodeIds.push_back(anAdvEdge.nodeId[0]);
        sNodeIds.push_back(anAdvEdge.nodeId[1]);

        std::sort(sNodeIds.begin(), sNodeIds.end());

        for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j) {

            // Discard same edge
            if (sEdges[j] == anAdvEdge.id)
                continue;

            SAdvancingFrontEdge& anotherEdge = m_anAdvFront[sEdges[j]];

            std::vector<Integer> sOtherNodeIds;

            sOtherNodeIds.push_back(anotherEdge.nodeId[0]);
            sOtherNodeIds.push_back(anotherEdge.nodeId[1]);

            std::sort(sOtherNodeIds.begin(), sOtherNodeIds.end());

            if (sOtherNodeIds[0] == sNodeIds[0] && sOtherNodeIds[1] == sNodeIds[1]) {

                aDuplicateEdgeId = sEdges[j];
                return true;
            }
        }

        return false;
    }

    template <typename Real>
    bool CMshTriangleMesher<Real>::edgeOk(SAdvancingFrontEdge& anAdvEdge, CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, std::vector<Integer>& sEdges, const Real aTolerance)
    {

        CGeoCoordinate<Real> aPoint;

        CGeoLine<Real> aLine1(aNode1, aNode2);

        for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j) {

            // Exclude current edge
            if (sEdges[j] == anAdvEdge.id)
                continue;

            SAdvancingFrontEdge& anotherEdge = m_anAdvFront[sEdges[j]];

            if (anotherEdge.remove)
                continue;

            CGeoLine<Real>& aLine2 = anotherEdge.line;

            CGeoIntersectionType anIntersectionType;

            if (aLine1.intersects(aLine2, aPoint, anIntersectionType, aTolerance)) {
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

        for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j) {

            if (m_anAdvFront[sEdges[j]].remove)
                continue;

            for (Integer k = 0; k < 2; ++k) {

                Integer aNodeId = m_anAdvFront[sEdges[j]].nodeId[k];

                if (std::find(sNodes.begin(), sNodes.end(), aNodeId) == sNodes.end())
                    sNodes.push_back(aNodeId);
            }
        }

        for (Integer j = 0; j < static_cast<Integer>(m_interiorNodes.size()); ++j) {

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

        for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j) {

            if (sEdges[j] == anAdvEdgeId)
                continue;

            if (m_anAdvFront[sEdges[j]].remove)
                continue;

            CGeoCoordinate<Real> aNewPoint;

            CGeoIntersectionType anIntersectionType;

            if (m_anAdvFront[sEdges[j]].line.intersects(aLine, aNewPoint, anIntersectionType, aTolerance)) {

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

        for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j) {

            aNodeId = sNodes[j];

            CMshNode<Real>& aNode = m_surfaceMesh.node(aNodeId);

            if ((aNode - aNode1).norm() < aTolerance || (aNode - aNode2).norm() < aTolerance || (aNode - aNode3).norm() < aTolerance)
                continue;

            CGeoIntersectionType anIntersectionType;

            if (aTriangle.contains(aNode, anIntersectionType, aTolerance)) {

                if (anIntersectionType == IT_INTERNAL || anIntersectionType == IT_EDGE)
                    return true;
            }
        }

        return false;
    }

    template <typename Real>
    bool CMshTriangleMesher<Real>::checkDelaunay(CMshNode<Real>& aNewNode, const Real aTolerance)
    {

        for (Integer i = 0; i < m_surfaceMesh.nbElements(); ++i) {

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

        try {

            for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {

                Integer anAdvEdgeId = i;

                SAdvancingFrontEdge& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                if (anAdvEdge.remove)
                    continue;

                if (anAdvEdge.nodeId[0] == anAdvEdge.nodeId[1]) {
                    anAdvEdge.remove = true;
                    continue;
                }

                if (anAdvEdge.neighborId[0] >= static_cast<Integer>(m_anAdvFront.size()) || anAdvEdge.neighborId[1] >= static_cast<Integer>(m_anAdvFront.size())) {
                    throw std::out_of_range("Connectivity is out of range!");
                }

                SAdvancingFrontEdge& aPrevEdge = m_anAdvFront[anAdvEdge.neighborId[0]];
                SAdvancingFrontEdge& aNextEdge = m_anAdvFront[anAdvEdge.neighborId[1]];

                if (aPrevEdge.remove || aNextEdge.remove) {

                    for (Integer j = 0; j < static_cast<Integer>(sEdges.size()); ++j) {

                        Integer anotherEdgeId = sEdges[j];

                        // Discard same edge
                        if (anotherEdgeId == anAdvEdge.id)
                            continue;

                        SAdvancingFrontEdge& anotherEdge = m_anAdvFront[anotherEdgeId];

                        if (anotherEdge.remove)
                            continue;

                        if (anAdvEdge.nodeId[0] == anotherEdge.nodeId[1]) {
                            anAdvEdge.neighborId[0] = anotherEdgeId;
                            anotherEdge.neighborId[1] = anAdvEdgeId;
                        } else if (anAdvEdge.nodeId[1] == anotherEdge.nodeId[0]) {
                            anAdvEdge.neighborId[1] = anotherEdgeId;
                            anotherEdge.neighborId[0] = anAdvEdgeId;
                        }
                    }
                }
            }

        } catch (const std::exception& e) {
            std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        } catch (...) {
            std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        }
    }

    template <typename Real>
    void CMshTriangleMesher<Real>::cleanDuplicateEdges(std::vector<Integer>& sEdges, const Real aTolerance)
    {

        try {

            if (m_nextEdgeId < 2)
                return;

            bool bAdjustConnectivity = false;

            // Check last two edges
            for (Integer i = 0; i < 2; ++i) {

                Integer anAdvEdgeId = m_nextEdgeId - i - 1;

                SAdvancingFrontEdge& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                if (anAdvEdge.remove)
                    continue;

                Integer aDuplicateEdgeId;

                if (edgeExists(anAdvEdge, aDuplicateEdgeId, sEdges)) {

                    SAdvancingFrontEdge& aDuplicateEdge = m_anAdvFront[aDuplicateEdgeId];

                    this->removeEdge(anAdvEdge, aTolerance);
                    this->removeEdge(aDuplicateEdge, aTolerance);

                    bAdjustConnectivity = true;
                }
            }

            if (bAdjustConnectivity)
                this->adjustConnectivity(sEdges);

        } catch (const std::exception& e) {
            std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        } catch (...) {
            std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        }
    }

    template <typename Real>
    void CMshTriangleMesher<Real>::addTriangle(SAdvancingFrontEdge& anAdvEdge, const Integer aNodeId, std::vector<Integer>& sEdges, const Real aTolerance)
    {

        try {

            Integer aNodeId1 = anAdvEdge.nodeId[0];
            Integer aNodeId2 = anAdvEdge.nodeId[1];
            Integer aNodeId3 = aNodeId;

            if (aNodeId3 == aNodeId1 || aNodeId3 == aNodeId2) {
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
            aNewEdge1.triangleId = aNewElementId;
            aNewEdge1.nodeNotId3 = aNodeId2;

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
            aNewEdge2.triangleId = aNewElementId;
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

        } catch (const std::exception& e) {
            std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        } catch (...) {
            std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        }
    }

    template <typename Real>
    bool CMshTriangleMesher<Real>::remesh(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, Real meshSize)
    {

        std::stringstream ss(std::stringstream::in | std::stringstream::out);

        ss << meshSize;

        ENigMA::analytical::CAnaFunction<Real> aAnaFunction;

        aAnaFunction.set(ss.str());

        return this->remesh(anEdgeMesh, aAnaFunction);
    }

    template <typename Real>
    bool CMshTriangleMesher<Real>::remesh(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc)
    {

        Real x, y;

        meshSizeFunc.removeAllVariables();

        meshSizeFunc.defineVariable("x", x);
        meshSizeFunc.defineVariable("y", y);

        // Split edge mesh according to local mesh size
        Integer nbEdges = anEdgeMesh.nbElements();

        for (Integer i = 0; i < nbEdges; ++i) {

            Integer anAdvEdgeId = anEdgeMesh.elementId(i);

            CMshElement<Real>& anElement = anEdgeMesh.element(anAdvEdgeId);

            Integer aNodeId1 = anElement.nodeId(0);
            Integer aNodeId2 = anElement.nodeId(1);

            CMshNode<Real> aNode1 = anEdgeMesh.node(aNodeId1);
            CMshNode<Real> aNode2 = anEdgeMesh.node(aNodeId2);

            bool bConnectivity = true;

            Integer aFaceId1 = 0;
            Integer aFaceId2 = 0;

            if (anElement.nbFaceIds() == 2) {
                aFaceId1 = anElement.faceId(0);
                aFaceId2 = anElement.faceId(1);
            } else
                bConnectivity = false;

            CGeoVector<Real> v = aNode2 - aNode1;

            CMshNode<Real> aMidNode = (aNode1 + aNode2) * 0.5;

            x = aMidNode.x();
            y = aMidNode.y();

            Real localMeshSize = meshSizeFunc.evaluate();

            Integer ne = static_cast<Integer>(floor(v.norm() / localMeshSize + 0.5));

            if (ne > 0) {

                Real de = v.norm() / ne;

                CGeoNormal<Real> n = v;

                n.normalize();

                Integer aPrevNodeId = aNodeId1;
                Integer aPrevFaceId = aFaceId1;

                for (Integer j = 0; j < ne - 1; ++j) {

                    // Create new node
                    CMshNode<Real> aNode = aNode1 + (j + 1) * de * n;

                    Integer aNewNodeId = anEdgeMesh.nextNodeId();
                    anEdgeMesh.addNode(aNewNodeId, aNode);

                    CMshElement<Real> aNewElement(ET_BEAM);
                    aNewElement.addNodeId(aPrevNodeId);
                    aNewElement.addNodeId(aNewNodeId);

                    Integer aNewElementId = anEdgeMesh.nextElementId();

                    if (bConnectivity) {

                        std::vector<CMshFace<Real>> sFaces;
                        aNewElement.generateFaces(sFaces);

                        for (Integer j = 0; j < static_cast<Integer>(sFaces.size()); ++j) {

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

                    if (bConnectivity) {
                        anEdgeMesh.element(anAdvEdgeId).setFaceId(0, aPrevFaceId);
                        anEdgeMesh.element(anAdvEdgeId).setFaceId(1, aFaceId2);
                    }
                }
            }
        }

        Integer aFirstFaceId = 0;

        for (Integer i = 0; i < anEdgeMesh.nbFaces(); ++i) {

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
    bool CMshTriangleMesher<Real>::generate(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, Real meshSize, Real minQuality, const Real aTolerance)
    {

        try {

            std::vector<CGeoCoordinate<double>> sInteriorPoints;

            ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
            aAnaFunction.set(meshSize);

            return this->generate(anEdgeMesh, maxNbElements, sInteriorPoints, aAnaFunction, minQuality, aTolerance);

        } catch (const std::exception& e) {
            std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        } catch (...) {
            std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        }
    }

    template <typename Real>
    bool CMshTriangleMesher<Real>::generate(CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, Real meshSize, Real minQuality, const Real aTolerance)
    {

        try {

            ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
            aAnaFunction.set(meshSize);

            return this->generate(anEdgeMesh, maxNbElements, sInteriorPoints, aAnaFunction, minQuality, aTolerance);

        } catch (const std::exception& e) {
            std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            return false;
        } catch (...) {
            std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            return false;
        }
    }

    template <typename Real>
    bool CMshTriangleMesher<Real>::generate(CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minQuality, const Real aTolerance)
    {

        try {

            m_previousNbElements = 0;

            m_bStop = false;

            // Add boundary nodes to surface mesh
            m_surfaceMesh.reset();
            m_boundingBox.reset();

            for (Integer i = 0; i < anEdgeMesh.nbNodes(); ++i) {

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

            m_nextEdgeId = 0;

            for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i) {

                Integer anAdvEdgeId = anEdgeMesh.elementId(i);

                CMshElement<Real>& anElement = anEdgeMesh.element(anAdvEdgeId);

                if (anElement.elementType() == ET_BEAM)
                    newEdgeIds[anAdvEdgeId] = m_nextEdgeId++;
            }

            m_nextEdgeId = 0;

            for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i) {

                Integer anAdvEdgeId = anEdgeMesh.elementId(i);

                CMshElement<Real>& anElement = anEdgeMesh.element(anAdvEdgeId);

                if (anElement.elementType() == ET_BEAM) {

                    SAdvancingFrontEdge anAdvEdge;

                    anAdvEdge.id = m_nextEdgeId++;
                    anAdvEdge.remove = false;
                    anAdvEdge.boundary = true;

                    for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                        anAdvEdge.nodeId[j] = anElement.nodeId(j);

                    anAdvEdge.triangleId = std::numeric_limits<Integer>::max();
                    anAdvEdge.nodeNotId3 = std::numeric_limits<Integer>::max();

                    anAdvEdge.build(anEdgeMesh);

                    for (Integer j = 0; j < anElement.nbFaceIds(); ++j) {

                        anAdvEdge.neighborId[j] = std::numeric_limits<Integer>::max();

                        Integer aFaceId = anElement.faceId(j);
                        CMshFace<Real>& aFace = anEdgeMesh.face(aFaceId);

                        if (aFace.hasPair()) {
                            Integer aPairFaceId = aFace.pairFaceId();
                            Integer anElementId = anEdgeMesh.face(aPairFaceId).elementId();

                            if (newEdgeIds.find(anElementId) != newEdgeIds.end()) {
                                anAdvEdge.neighborId[j] = newEdgeIds.at(anElementId);
                            } else
                                std::cout << "Error: element id = " << anElementId << " not found!" << std::endl;

                        } else {
                            throw std::runtime_error("Boundary is open!");
                        }
                    }

                    m_anAdvFront.push_back(anAdvEdge);

                    this->addEdgeToRtree(anAdvEdge, aTolerance);
                }
            }

            // Add interior nodes
            for (Integer i = 0; i < static_cast<Integer>(sInteriorPoints.size()); ++i) {

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

            res ? res = this->advancingFrontTriMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 0.75, 0.10, false, 0, aTolerance) : res = false;
            res ? res = this->advancingFrontTriMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 1.50, 0.00, false, 0, aTolerance) : res = false;

            m_surfaceMesh.removeDanglingNodes();
            m_surfaceMesh.renumber();

            if (this->frontSize() == 0) {
                m_surfaceMesh.generateFaces(aTolerance);
                return true;
            } else {
                std::cout << "Meshing error!" << std::endl;
                return false;
            }

        } catch (const std::exception& e) {
            std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            return false;
        } catch (...) {
            std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            return false;
        }
    }

    template <typename Real>
    Integer CMshTriangleMesher<Real>::frontSize()
    {

        Integer n = 0;

        for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {
            if (!m_anAdvFront[i].remove)
                n++;
        }

        return n;
    }

    template <typename Real>
    bool CMshTriangleMesher<Real>::advancingFrontTriMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real sizeFactor, Real shrinkFactor, Real expandFactor, Real minQuality, const bool bCheckDelaunay, Integer firstIndex, const Real aTolerance)
    {

        bool res = true;

        try {

            Real x, y;

            meshSizeFunc.removeAllVariables();

            meshSizeFunc.defineVariable("x", x);
            meshSizeFunc.defineVariable("y", y);

            const Real pi = std::acos(-1.0);

            Real sumMeshSize = 0;
            Real averageMeshSize = 0;
            Integer nMeshSize = 0;

            for (Integer i = firstIndex; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {

                this->checkUpdate();

                if (m_bStop)
                    return false;

                if (maxNbElements > 0) {
                    if (m_surfaceMesh.nbElements() >= maxNbElements) {
                        std::cout << "Max number of elements (" << maxNbElements << ") reached!" << std::endl;
                        return false;
                    }
                }

                if (!m_anAdvFront[i].remove) {

                    Integer anAdvEdgeId = i;

                    SAdvancingFrontEdge& anAdvEdge = m_anAdvFront[anAdvEdgeId];

                    Integer aNodeId1 = anAdvEdge.nodeId[0];
                    Integer aNodeId2 = anAdvEdge.nodeId[1];

                    CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
                    CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);

                    Integer aNodeId3 = m_anAdvFront[anAdvEdge.neighborId[0]].nodeId[0];
                    Integer aNodeId4 = m_anAdvFront[anAdvEdge.neighborId[1]].nodeId[1];

                    CMshNode<Real> aMidNode = (aNode1 + aNode2) * 0.5;

                    x = aMidNode.x();
                    y = aMidNode.y();

                    // Rotate vector by 90�
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

                    v.rotate(pi * 0.5);

                    // Equilateral triangle (height to edge ratio)
                    v *= sqrt(3.0) * 0.5;

                    // Add point to form triangle with correct spacing
                    CMshNode<Real> aNewNode = aMidNode + v * localMeshSize * sizeFactor;
                    Integer aNewNodeId = m_surfaceMesh.nextNodeId();

                    // Get closest edges
                    CGeoBoundingBox<Real> aBoundingBox;
                    aBoundingBox.addCoordinate(aNode1);
                    aBoundingBox.addCoordinate(aNode2);
                    aBoundingBox.addCoordinate(aNewNode);
                    aBoundingBox.grow(localMeshSize * sizeFactor * 0.5);

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

                    if (aNodeId3 == aNodeId4) {

                        CMshNode<Real>& aNode3 = m_surfaceMesh.node(aNodeId3);

                        CMshTriangle<Real> aNewTriangle;

                        aNewTriangle.addVertex(aNode1);
                        aNewTriangle.addVertex(aNode2);
                        aNewTriangle.addVertex(aNode3);

                        aNewTriangle.calculateArea();

                        if (aNewTriangle.normal().z() > aTolerance) {

                            Integer wNodeId;

                            if (this->edgeOk(anAdvEdge, aNode1, aNode3, sEdges, aTolerance) && this->edgeOk(anAdvEdge, aNode2, aNode3, sEdges, aTolerance) && !this->triangleContainsNode(aNode1, aNode2, aNode3, wNodeId, sNodes, aTolerance)) {

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

                    for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j) {

                        Integer aNodeId = sNodes[j];

                        if (aNodeId == aNodeId1 || aNodeId == aNodeId2)
                            continue;

                        CMshNode<Real>& aNode = m_surfaceMesh.node(aNodeId);

                        Real d = (aNode - aNewNode).norm();

                        Real factor = 1.0;

                        if (aNodeId == aNodeId3) {
                            Real angle = anAdvEdge.line.vector().angle(m_anAdvFront[anAdvEdge.neighborId[0]].line.vector());

                            if (angle > pi * 0.5)
                                factor = 1.5;
                        }

                        if (aNodeId == aNodeId4) {
                            Real angle = anAdvEdge.line.vector().angle(m_anAdvFront[anAdvEdge.neighborId[1]].line.vector());

                            if (angle > pi * 0.5)
                                factor = 1.5;
                        }

                        // Use closest node
                        if (d < localMeshSize * sizeFactor * expandFactor * factor) {

                            CMshTriangle<Real> aNewTriangle1;

                            aNewTriangle1.addVertex(aNode1);
                            aNewTriangle1.addVertex(aNode2);
                            aNewTriangle1.addVertex(aNode);

                            aNewTriangle1.calculateArea();

                            if (aNewTriangle1.normal().z() > aTolerance) {

                                aNewTriangle1.calculateQuality();

                                Real q1 = aNewTriangle1.quality();

                                if (q1 > qmax)
                                    q1 += this->edgeOk(anAdvEdge, aNode1, aNode, sEdges, aTolerance) ? 0.0 : -2.0;

                                if (q1 > qmax)
                                    q1 += this->edgeOk(anAdvEdge, aNode2, aNode, sEdges, aTolerance) ? 0.0 : -2.0;

                                if (q1 > qmax) {
                                    Integer wNodeId;
                                    q1 += !this->triangleContainsNode(aNode1, aNode2, aNode, wNodeId, sNodes, aTolerance) ? 0.0 : -2.0;
                                }

                                if (q1 > qmax) {
                                    qmax = q1;
                                    anExistingNodeId = aNodeId;
                                    anExistingNode = aNode;
                                }
                            }
                        }
                    }

                    if (qmax > minQuality) {

                        for (Integer k = 0; k < static_cast<Integer>(m_interiorNodes.size()); ++k) {

                            if (m_interiorNodes[k].remove)
                                continue;

                            if (m_interiorNodes[k].nodeId == anExistingNodeId)
                                m_interiorNodes[k].remove = true;
                        }

                        this->addTriangle(anAdvEdge, anExistingNodeId, sEdges, aTolerance);
                        res = true;

                    } else {

                        CGeoLine<Real> aLine(aMidNode, aNewNode);

                        Real dmin = findShortestDistance(sEdges, aLine, anAdvEdgeId, aTolerance);

                        if (dmin > localMeshSize * sizeFactor * shrinkFactor * 0.25 && dmin < localMeshSize * sizeFactor * expandFactor)
                            aNewNode = aMidNode + v * dmin * 0.5;

                        CMshTriangle<Real> aNewTriangle2;

                        aNewTriangle2.addVertex(aNode1);
                        aNewTriangle2.addVertex(aNode2);
                        aNewTriangle2.addVertex(aNewNode);

                        aNewTriangle2.calculateArea();

                        if (aNewTriangle2.normal().z() > aTolerance) {

                            aNewTriangle2.calculateQuality();

                            Real q2 = aNewTriangle2.quality();

                            if (q2 > 1.5 * minQuality)
                                q2 += this->edgeOk(anAdvEdge, aNode1, aNewNode, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q2 > 1.5 * minQuality)
                                q2 += this->edgeOk(anAdvEdge, aNode2, aNewNode, sEdges, aTolerance) ? 0.0 : -2.0;

                            if (q2 > 1.5 * minQuality) {
                                Integer wNodeId;
                                q2 += !this->triangleContainsNode(aNode1, aNode2, aNewNode, wNodeId, sNodes, aTolerance) ? 0.0 : -2.0;
                            }

                            if (q2 > 1.5 * minQuality && bCheckDelaunay)
                                q2 += this->checkDelaunay(aNewNode, aTolerance) ? 0.0 : -2.0;

                            if (q2 > 1.5 * minQuality) {

                                // Create a new node
                                m_surfaceMesh.addNode(aNewNodeId, aNewNode);

                                this->addTriangle(anAdvEdge, aNewNodeId, sEdges, aTolerance);
                                res = true;

                                if (!m_boundingBox.contains(aNewNode, aTolerance)) {
                                    throw std::runtime_error("Node is outside boundary!");
                                }
                            }
                        }
                    }
                }
            }

            if (onUpdate != nullptr)
                onUpdate(0);

        } catch (const std::exception& e) {
            std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        } catch (...) {
            std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            throw;
        }

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
    void CMshTriangleMesher<Real>::applyFixedBoundary(CMshMesh<Real>& anEdgeMesh, const Real aTolerance)
    {

        // Discover double faces
        CGeoHashGrid<Real> aHashGrid;

        std::vector<CGeoCoordinate<Real>> sCenterCoordinates;

        for (Integer i = 0; i < m_surfaceMesh.nbFaces(); ++i) {

            Integer aFaceId = m_surfaceMesh.faceId(i);
            CMshFace<Real>& aFace = m_surfaceMesh.face(aFaceId);

            CGeoCoordinate<Real> aCenterCoordinate(0.0, 0.0, 0.0);

            for (Integer j = 0; j < aFace.nbNodeIds(); ++j) {

                Integer aNodeId = aFace.nodeId(j);
                CMshNode<Real>& aNode = m_surfaceMesh.node(aNodeId);

                aCenterCoordinate += aNode;
            }

            if (aFace.nbNodeIds() > 0)
                aCenterCoordinate /= static_cast<Real>(aFace.nbNodeIds());

            sCenterCoordinates.push_back(aCenterCoordinate);

            aHashGrid.addGeometricObject(aFaceId, aCenterCoordinate);
        }

        aHashGrid.build();

        for (Integer i = 0; i < anEdgeMesh.nbElements(); ++i) {

            Integer anElementId = anEdgeMesh.elementId(i);

            CMshElement<Real>& anElement = anEdgeMesh.element(anElementId);

            CGeoCoordinate<Real> aCenterCoordinate(0.0, 0.0, 0.0);

            for (Integer j = 0; j < anElement.nbNodeIds(); ++j) {

                Integer aNodeId = anElement.nodeId(j);
                CMshNode<Real> aNode = anEdgeMesh.node(aNodeId);

                aCenterCoordinate += aNode;
            }

            if (anElement.nbNodeIds() > 0)
                aCenterCoordinate /= static_cast<Real>(anElement.nbNodeIds());

            std::vector<Integer> sCoordinates;

            aHashGrid.find(sCoordinates, aCenterCoordinate, aTolerance);

            for (Integer j = 0; j < static_cast<Integer>(sCoordinates.size()); ++j) {

                Integer aFaceId = sCoordinates[j];
                m_surfaceMesh.face(aFaceId).setHasPair(false);
            }
        }
    }

    template <typename Real>
    void CMshTriangleMesher<Real>::flipEdges(const Real aTolerance)
    {

        std::map<Integer, bool> sFlipped;

        for (Integer i = 0; i < m_surfaceMesh.nbElements(); ++i) {

            Integer anElementId = m_surfaceMesh.elementId(i);
            sFlipped[anElementId] = false;
        }

        for (Integer i = 0; i < m_surfaceMesh.nbFaces(); ++i) {

            Integer aFaceId = m_surfaceMesh.faceId(i);
            CMshFace<Real> aFace = m_surfaceMesh.face(aFaceId);

            if (aFace.faceType() != FT_LINE)
                continue;

            Integer anElementId = aFace.elementId();

            if (sFlipped.at(anElementId))
                continue;

            CMshElement<Real>& anElement = m_surfaceMesh.element(anElementId);

            if (anElement.elementType() != ET_TRIANGLE)
                continue;

            Integer aNodeId1 = aFace.nodeId(0);
            Integer aNodeId2 = aFace.nodeId(1);

            Integer aNodeId3;

            bool bNodeId3 = false;

            for (Integer j = 0; j < anElement.nbNodeIds(); ++j) {

                if (anElement.nodeId(j) != aNodeId1 && anElement.nodeId(j) != aNodeId2) {
                    aNodeId3 = anElement.nodeId(j);
                    bNodeId3 = true;
                    break;
                }
            }

            if (!bNodeId3)
                continue;

            if (aFace.hasPair()) {

                Integer pairFaceId = aFace.pairFaceId();
                CMshFace<Real> aPairFace = m_surfaceMesh.face(pairFaceId);

                Integer aNeighborId = aPairFace.elementId();

                if (sFlipped.at(aNeighborId))
                    continue;

                CMshElement<Real>& aNeighbor = m_surfaceMesh.element(aNeighborId);

                if (aNeighbor.elementType() != ET_TRIANGLE)
                    continue;

                Integer aNodeId4;

                bool bNodeId4 = false;

                for (Integer j = 0; j < aNeighbor.nbNodeIds(); ++j) {

                    if (aNeighbor.nodeId(j) != aNodeId1 && aNeighbor.nodeId(j) != aNodeId2) {
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

                if (std::min(q3, q4) > std::min(q1, q2) && aTriangle3.normal().z() > aTolerance && aTriangle4.normal().z() > aTolerance) {

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
    void CMshTriangleMesher<Real>::relaxNodes(const Real aTolerance)
    {

        std::map<Integer, bool> bBoundaryNode;

        for (Integer i = 0; i < m_surfaceMesh.nbNodes(); ++i) {

            Integer aNodeId = m_surfaceMesh.nodeId(i);
            bBoundaryNode[aNodeId] = false;
        }

        for (Integer i = 0; i < m_surfaceMesh.nbFaces(); ++i) {

            Integer aFaceId = m_surfaceMesh.faceId(i);
            CMshFace<Real>& aFace = m_surfaceMesh.face(aFaceId);

            if (aFace.hasPair())
                continue;

            for (Integer j = 0; j < aFace.nbNodeIds(); ++j) {

                Integer aNodeId = aFace.nodeId(j);
                bBoundaryNode[aNodeId] = true;
            }
        }

        CMshMeshQuery<Real> aMeshQuery(m_surfaceMesh);

        std::vector<Integer> sElementIds;

        for (Integer i = 0; i < m_surfaceMesh.nbNodes(); ++i) {

            Integer aMovingNodeId = m_surfaceMesh.nodeId(i);
            CMshNode<Real>& aMovingNode = m_surfaceMesh.node(aMovingNodeId);

            if (bBoundaryNode.at(aMovingNodeId))
                continue;

            aMeshQuery.elementsSharingNode(aMovingNodeId, sElementIds);

            std::vector<Integer> sNodeIds;
            sNodeIds.push_back(aMovingNodeId);

            for (Integer j = 0; j < static_cast<Integer>(sElementIds.size()); ++j) {

                Integer anElementId = sElementIds[j];

                CMshElement<Real>& anElement = m_surfaceMesh.element(anElementId);

                if (anElement.elementType() != ET_TRIANGLE)
                    continue;

                for (Integer k = 0; k < anElement.nbNodeIds(); ++k) {

                    Integer aNodeId = anElement.nodeId(k);

                    if (std::find(sNodeIds.begin(), sNodeIds.end(), aNodeId) == sNodeIds.end())
                        sNodeIds.push_back(aNodeId);
                }
            }

            if (sNodeIds.size() > 1) {

                CMshNode<Real> aNewNode(0.0, 0.0, 0.0);

                for (Integer k = 0; k < static_cast<Integer>(sNodeIds.size()); ++k)
                    aNewNode += m_surfaceMesh.node(sNodeIds[k]);

                aNewNode /= static_cast<Real>(sNodeIds.size());

                Real sumq1 = 0.0;
                Real sumq2 = 0.0;

                bool bMoveNode = true;

                for (Integer k = 0; k < static_cast<Integer>(sElementIds.size()); ++k) {

                    CMshElement<Real>& aModElement = m_surfaceMesh.element(sElementIds[k]);

                    if (aModElement.elementType() != ET_TRIANGLE)
                        continue;

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

                    if (aTriangle2.normal().z() <= aTolerance) {
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
    void CMshTriangleMesher<Real>::collapseEdges(Real collapseSize, const Real aTolerance)
    {

        std::map<Integer, bool> bBoundaryNode;

        for (Integer i = 0; i < m_surfaceMesh.nbNodes(); ++i) {

            Integer aNodeId = m_surfaceMesh.nodeId(i);
            bBoundaryNode[aNodeId] = false;
        }

        for (Integer i = 0; i < m_surfaceMesh.nbFaces(); ++i) {

            Integer aFaceId = m_surfaceMesh.faceId(i);
            CMshFace<Real>& aFace = m_surfaceMesh.face(aFaceId);

            if (aFace.hasPair())
                continue;

            for (Integer j = 0; j < aFace.nbNodeIds(); ++j) {

                Integer aNodeId = aFace.nodeId(j);
                bBoundaryNode[aNodeId] = true;
            }
        }

        for (Integer i = 0; i < m_surfaceMesh.nbElements(); ++i) {

            Integer anElementId = m_surfaceMesh.elementId(i);
            CMshElement<Real>& anElement = m_surfaceMesh.element(anElementId);

            if (anElement.elementType() != ET_TRIANGLE)
                continue;

            for (Integer j = 0; j < anElement.nbNodeIds(); ++j) {

                Integer aNodeId1 = anElement.nodeId((j + 0) % 3);
                Integer aNodeId2 = anElement.nodeId((j + 1) % 3);

                CMshNode<Real>& aNode1 = m_surfaceMesh.node(aNodeId1);
                CMshNode<Real>& aNode2 = m_surfaceMesh.node(aNodeId2);

                Real anEdgeLength = (aNode1 - aNode2).norm();

                if (bBoundaryNode.at(aNodeId2))
                    continue;

                if (anEdgeLength < collapseSize * 0.75) {
                    // Try to collapse
                    aNode2 = aNode1;
                    anElement.setElementType(ET_NONE);
                    break;
                }
            }
        }

        // Delete invalid elements
        for (Integer i = 0; i < m_surfaceMesh.nbElements(); ++i) {

            Integer anElementId = m_surfaceMesh.elementId(i);
            CMshElement<Real>& anElement = m_surfaceMesh.element(anElementId);

            if (anElement.elementType() == ET_NONE) {
                m_surfaceMesh.removeElement(anElementId);
            }
        }

        m_surfaceMesh.mergeNodes(aTolerance);
        m_surfaceMesh.renumber();

        for (Integer i = 0; i < m_surfaceMesh.nbElements(); ++i) {

            Integer anElementId = m_surfaceMesh.elementId(i);
            CMshElement<Real>& anElement = m_surfaceMesh.element(anElementId);

            if (anElement.elementType() == ET_NONE) {
                m_surfaceMesh.element(anElementId).setElementType(ET_TRIANGLE);
            }
        }

        m_surfaceMesh.generateFaces(aTolerance);
    }
}
}
