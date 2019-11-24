// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoNormal.hpp"
#include "GeoSphere.hpp"
#include "MshMeshQuery.hpp"
#include "MshTetrahedron.hpp"

using namespace ENigMA::geometry;

namespace ENigMA {

namespace mesh {

    template <typename Real>
    CMshTetrahedronMesher<Real>::CMshTetrahedronMesher()
        : m_timeInterval(1)
        , m_dataInterval(100)
        , m_begin(0)
        , m_end(0)
        , m_bStop(false)
        , onUpdate(nullptr)
        , m_previousNbElements(0)
    {
    }

    template <typename Real>
    CMshTetrahedronMesher<Real>::~CMshTetrahedronMesher()
    {
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::checkUpdate()
    {

        if (onUpdate != nullptr) {
            m_end = clock();

            if (static_cast<Real>(m_end - m_begin) / CLOCKS_PER_SEC > m_timeInterval) {
                m_begin = m_end;

                if (m_volumeMesh.nbElements() > m_previousNbElements + m_dataInterval) {
                    m_previousNbElements = m_volumeMesh.nbElements();
                    onUpdate(true);
                } else
                    onUpdate(false);
            }
        }
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::removeTriangle(SAdvancingFrontTriangle& anAdvTriangle, const Real aTolerance)
    {

        anAdvTriangle.remove = true;
        this->removeTriangleFromRtree(anAdvTriangle, aTolerance);
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::addTriangleToRtree(SAdvancingFrontTriangle& anAdvTriangle, const Real aTolerance)
    {

        Integer aNodeId1 = anAdvTriangle.nodeId[0];
        Integer aNodeId2 = anAdvTriangle.nodeId[1];
        Integer aNodeId3 = anAdvTriangle.nodeId[2];

        CMshNode<Real>& aNode1 = m_volumeMesh.node(aNodeId1);
        CMshNode<Real>& aNode2 = m_volumeMesh.node(aNodeId2);
        CMshNode<Real>& aNode3 = m_volumeMesh.node(aNodeId3);

        CGeoBoundingBox<Real> aBoundingBox;

        aBoundingBox.addCoordinate(aNode1);
        aBoundingBox.addCoordinate(aNode2);
        aBoundingBox.addCoordinate(aNode3);

        aBoundingBox.grow(aTolerance);

        m_tree.addGeometricObject(anAdvTriangle.id, aBoundingBox);
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::removeTriangleFromRtree(SAdvancingFrontTriangle& anAdvTriangle, const Real aTolerance)
    {

        Integer aNodeId1 = anAdvTriangle.nodeId[0];
        Integer aNodeId2 = anAdvTriangle.nodeId[1];
        Integer aNodeId3 = anAdvTriangle.nodeId[2];

        CGeoBoundingBox<Real> aBoundingBox;

        aBoundingBox.addCoordinate(m_volumeMesh.node(aNodeId1));
        aBoundingBox.addCoordinate(m_volumeMesh.node(aNodeId2));
        aBoundingBox.addCoordinate(m_volumeMesh.node(aNodeId3));

        aBoundingBox.grow(aTolerance);

        m_tree.removeGeometricObject(anAdvTriangle.id, aBoundingBox);
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::triangleExists(SAdvancingFrontTriangle& anAdvTriangle, Integer& aDuplicateTriangleId, std::vector<Integer>& sTriangles)
    {

        std::vector<Integer> sNodeIds;

        sNodeIds.push_back(anAdvTriangle.nodeId[0]);
        sNodeIds.push_back(anAdvTriangle.nodeId[1]);
        sNodeIds.push_back(anAdvTriangle.nodeId[2]);

        std::sort(sNodeIds.begin(), sNodeIds.end());

        for (Integer j = 0; j < static_cast<Integer>(sTriangles.size()); ++j) {

            // Discard same triangle
            if (sTriangles[j] == anAdvTriangle.id)
                continue;

            SAdvancingFrontTriangle& anotherTriangle = m_anAdvFront[sTriangles[j]];

            std::vector<Integer> sOtherNodeIds;

            sOtherNodeIds.push_back(anotherTriangle.nodeId[0]);
            sOtherNodeIds.push_back(anotherTriangle.nodeId[1]);
            sOtherNodeIds.push_back(anotherTriangle.nodeId[2]);

            std::sort(sOtherNodeIds.begin(), sOtherNodeIds.end());

            if (sOtherNodeIds[0] == sNodeIds[0] && sOtherNodeIds[1] == sNodeIds[1] && sOtherNodeIds[2] == sNodeIds[2]) {

                aDuplicateTriangleId = sTriangles[j];
                return true;
            }
        }

        return false;
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::triangleOk(SAdvancingFrontTriangle& anAdvTriangle, CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, std::vector<Integer>& sTriangles, const Real aTolerance)
    {

        CGeoTriangle<Real> aTriangle1;

        aTriangle1.addVertex(aNode1);
        aTriangle1.addVertex(aNode2);
        aTriangle1.addVertex(aNode3);

        for (Integer j = 0; j < static_cast<Integer>(sTriangles.size()); ++j) {

            // Exclude current triangle
            if (sTriangles[j] == anAdvTriangle.id)
                continue;

            SAdvancingFrontTriangle& anotherAdvTriangle = m_anAdvFront[sTriangles[j]];

            if (anotherAdvTriangle.remove)
                continue;

            CGeoTriangle<Real>& aTriangle2 = anotherAdvTriangle.triangle;

            CGeoIntersectionType anIntersectionType;

            if (aTriangle1.intersects(aTriangle2, anIntersectionType, aTolerance)) {

                if (anIntersectionType == IT_INTERNAL || anIntersectionType == IT_SWAP)
                    return false;
            }
        }

        return true;
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::findClosestNodes(std::vector<Integer>& sTriangles, std::vector<Integer>& sNodes)
    {

        for (Integer j = 0; j < static_cast<Integer>(sTriangles.size()); ++j) {

            if (m_anAdvFront[sTriangles[j]].remove)
                continue;

            for (Integer k = 0; k < 3; ++k) {

                Integer aNodeId = m_anAdvFront[sTriangles[j]].nodeId[k];

                if (std::find(sNodes.begin(), sNodes.end(), aNodeId) == sNodes.end())
                    sNodes.push_back(aNodeId);
            }
        }

        for (Integer j = 0; j < static_cast<Integer>(m_innerNodes.size()); ++j) {

            if (m_innerNodes[j].remove)
                continue;

            Integer aNodeId = m_innerNodes[j].nodeId;

            if (std::find(sNodes.begin(), sNodes.end(), aNodeId) == sNodes.end())
                sNodes.push_back(aNodeId);
        }
    }

    template <typename Real>
    Real CMshTetrahedronMesher<Real>::findShortestDistance(std::vector<Integer>& sTriangles, CGeoLine<Real>& aLine, Integer anAdvTriangleId, const Real aTolerance)
    {

        Real dmin = std::numeric_limits<Real>::max();

        for (Integer j = 0; j < static_cast<Integer>(sTriangles.size()); ++j) {

            if (sTriangles[j] == anAdvTriangleId)
                continue;

            if (m_anAdvFront[sTriangles[j]].remove)
                continue;

            CGeoCoordinate<Real> aNewPoint;

            CGeoIntersectionType anIntersectionType;

            if (m_anAdvFront[sTriangles[j]].triangle.intersects(aLine, aNewPoint, anIntersectionType, aTolerance)) {

                Real distance = (aNewPoint - aLine.startPoint()).norm();

                if (distance < dmin)
                    dmin = distance;
            }
        }

        return dmin;
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::tetrahedronContainsNode(CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, CMshNode<Real>& aNode4, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance)
    {

        CMshTetrahedron<Real> aTetrahedron;

        aTetrahedron.addVertex(aNode1);
        aTetrahedron.addVertex(aNode2);
        aTetrahedron.addVertex(aNode3);
        aTetrahedron.addVertex(aNode4);

        for (Integer j = 0; j < static_cast<Integer>(sNodes.size()); ++j) {

            aNodeId = sNodes[j];

            CMshNode<Real>& aNode = m_volumeMesh.node(aNodeId);

            if ((aNode - aNode1).norm() < aTolerance || (aNode - aNode2).norm() < aTolerance || (aNode - aNode3).norm() < aTolerance || (aNode - aNode4).norm() < aTolerance)
                continue;

            if (aTetrahedron.contains(aNode, aTolerance))
                return true;
        }

        return false;
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::checkDelaunay(CMshNode<Real>& aNewNode, const Real aTolerance)
    {

        for (Integer i = 0; i < m_volumeMesh.nbElements(); ++i) {

            Integer anElementId = m_volumeMesh.elementId(i);

            CMshElement<Real>& anElement = m_volumeMesh.element(anElementId);

            if (anElement.elementType() != ET_TETRAHEDRON)
                continue;

            Integer aNodeId1 = anElement.nodeId(0);
            Integer aNodeId2 = anElement.nodeId(1);
            Integer aNodeId3 = anElement.nodeId(2);
            Integer aNodeId4 = anElement.nodeId(3);

            CMshNode<Real>& aNode1 = m_volumeMesh.node(aNodeId1);
            CMshNode<Real>& aNode2 = m_volumeMesh.node(aNodeId2);
            CMshNode<Real>& aNode3 = m_volumeMesh.node(aNodeId3);
            CMshNode<Real>& aNode4 = m_volumeMesh.node(aNodeId4);

            CGeoSphere<Real> aSphere(aNode1, aNode2, aNode3, aNode4);

            if (aSphere.contains(aNewNode, aTolerance))
                return false;
        }

        return true;
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::pairEdges(SAdvancingFrontTriangle& anAdvTriangle1, SAdvancingFrontTriangle& anAdvTriangle2)
    {

        for (Integer k = 0; k < 3; ++k) {

            for (Integer l = 0; l < 3; l++) {

                if ((anAdvTriangle1.nodeId[(k + 0) % 3] == anAdvTriangle2.nodeId[(l + 0) % 3] && anAdvTriangle1.nodeId[(k + 1) % 3] == anAdvTriangle2.nodeId[(l + 1) % 3]) || (anAdvTriangle1.nodeId[(k + 0) % 3] == anAdvTriangle2.nodeId[(l + 1) % 3] && anAdvTriangle1.nodeId[(k + 1) % 3] == anAdvTriangle2.nodeId[(l + 0) % 3])) {
                    anAdvTriangle1.neighborId[k] = anAdvTriangle2.id;
                    anAdvTriangle1.nodeNotId[k] = (l + 2) % 3;

                    anAdvTriangle2.neighborId[l] = anAdvTriangle1.id;
                    anAdvTriangle2.nodeNotId[l] = (k + 2) % 3;
                    break;
                }
            }
        }

        return false;
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::adjustConnectivity(std::vector<Integer>& sTriangles)
    {
        for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {

            Integer anAdvTriangleId = i;

            SAdvancingFrontTriangle& anAdvTriangle = m_anAdvFront[anAdvTriangleId];

            if (anAdvTriangle.remove)
                continue;

            if (anAdvTriangle.nodeId[0] == anAdvTriangle.nodeId[1] || anAdvTriangle.nodeId[1] == anAdvTriangle.nodeId[2] || anAdvTriangle.nodeId[2] == anAdvTriangle.nodeId[0]) {
                std::cout << "Error: invalid front triangle!" << std::endl;
                anAdvTriangle.remove = true;
                continue;
            }

            if (anAdvTriangle.neighborId[0] >= static_cast<Integer>(m_anAdvFront.size()) || anAdvTriangle.neighborId[1] >= static_cast<Integer>(m_anAdvFront.size()) || anAdvTriangle.neighborId[2] >= static_cast<Integer>(m_anAdvFront.size())) {
                throw std::out_of_range("Connectivity is out of range!");
            }

            if (m_anAdvFront[anAdvTriangle.neighborId[0]].remove || m_anAdvFront[anAdvTriangle.neighborId[1]].remove || m_anAdvFront[anAdvTriangle.neighborId[2]].remove) {

                for (Integer j = 0; j < static_cast<Integer>(sTriangles.size()); ++j) {

                    Integer anotherAdvTriangleId = sTriangles[j];

                    // Discard same triangle
                    if (anotherAdvTriangleId == anAdvTriangle.id)
                        continue;

                    SAdvancingFrontTriangle& anotherAdvTriangle = m_anAdvFront[anotherAdvTriangleId];

                    if (anotherAdvTriangle.remove)
                        continue;

                    this->pairEdges(anAdvTriangle, anotherAdvTriangle);
                }

                if (anAdvTriangle.neighborId[0] == std::numeric_limits<Integer>::max() || anAdvTriangle.neighborId[1] == std::numeric_limits<Integer>::max() || anAdvTriangle.neighborId[2] == std::numeric_limits<Integer>::max()) {
                    throw std::runtime_error("Boundary is open!");
                }
            }
        }
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::cleanDuplicateTriangles(std::vector<Integer>& sTriangles, const Real aTolerance)
    {
        if (m_nextTriangleId < 3)
            return;

        bool bAdjustConnectivity = false;

        // Check last three triangles
        for (Integer i = 0; i < 3; ++i) {

            Integer anAdvTriangleId = m_nextTriangleId - i - 1;

            SAdvancingFrontTriangle& anAdvTriangle = m_anAdvFront[anAdvTriangleId];

            if (anAdvTriangle.remove)
                continue;

            Integer aDuplicateTriangleId;

            if (this->triangleExists(anAdvTriangle, aDuplicateTriangleId, sTriangles)) {

                SAdvancingFrontTriangle& aDuplicateTriangle = m_anAdvFront[aDuplicateTriangleId];

                this->removeTriangle(anAdvTriangle, aTolerance);
                this->removeTriangle(aDuplicateTriangle, aTolerance);

                bAdjustConnectivity = true;
            }
        }

        if (bAdjustConnectivity)
            this->adjustConnectivity(sTriangles);
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::addTetrahedron(SAdvancingFrontTriangle& anAdvTriangle, const Integer aNodeId, std::vector<Integer>& sTriangles, const Real aTolerance)
    {
        Integer aNodeId1 = anAdvTriangle.nodeId[0];
        Integer aNodeId2 = anAdvTriangle.nodeId[1];
        Integer aNodeId3 = anAdvTriangle.nodeId[2];
        Integer aNodeId4 = aNodeId;

        if (aNodeId4 == aNodeId1 || aNodeId4 == aNodeId2 || aNodeId4 == aNodeId3) {
            std::cout << "Error: invalid tetrahedron!" << std::endl;
            return;
        }

        CMshElement<Real> aNewElement(ET_TETRAHEDRON);
        aNewElement.addNodeId(aNodeId1);
        aNewElement.addNodeId(aNodeId2);
        aNewElement.addNodeId(aNodeId3);
        aNewElement.addNodeId(aNodeId4);

        Integer aNewElementId = m_volumeMesh.nextElementId();

        m_volumeMesh.addElement(aNewElementId, aNewElement);

        SAdvancingFrontTriangle& aTriangle1 = m_anAdvFront[anAdvTriangle.neighborId[0]];
        SAdvancingFrontTriangle& aTriangle2 = m_anAdvFront[anAdvTriangle.neighborId[1]];
        SAdvancingFrontTriangle& aTriangle3 = m_anAdvFront[anAdvTriangle.neighborId[2]];

        Integer aNewTriangleId1 = m_nextTriangleId++;
        Integer aNewTriangleId2 = m_nextTriangleId++;
        Integer aNewTriangleId3 = m_nextTriangleId++;

        // Add triangle 1
        SAdvancingFrontTriangle aNewTriangle1;
        aNewTriangle1.id = aNewTriangleId1;
        aNewTriangle1.remove = false;
        aNewTriangle1.boundary = false;
        aNewTriangle1.nodeId[0] = aNodeId1;
        aNewTriangle1.nodeId[1] = aNodeId2;
        aNewTriangle1.nodeId[2] = aNodeId4;
        aNewTriangle1.neighborId[0] = anAdvTriangle.neighborId[0];
        aNewTriangle1.neighborId[1] = aNewTriangleId2;
        aNewTriangle1.neighborId[2] = aNewTriangleId3;
        aNewTriangle1.nodeNotId[0] = anAdvTriangle.nodeNotId[0];
        aNewTriangle1.nodeNotId[1] = 1;
        aNewTriangle1.nodeNotId[2] = 0;
        aNewTriangle1.tetrahedronId = aNewElementId;
        aNewTriangle1.nodeNotId4 = aNodeId3;

        aNewTriangle1.build(m_volumeMesh);

        // Correct connectivity
        aTriangle1.neighborId[(anAdvTriangle.nodeNotId[0] + 1) % 3] = aNewTriangleId1;
        aTriangle1.nodeNotId[(anAdvTriangle.nodeNotId[0] + 1) % 3] = 2;

        // Add triangle to rtree
        this->addTriangleToRtree(aNewTriangle1, aTolerance);

        // Add triangle 2
        SAdvancingFrontTriangle aNewTriangle2;
        aNewTriangle2.id = aNewTriangleId2;
        aNewTriangle2.remove = false;
        aNewTriangle2.boundary = false;
        aNewTriangle2.nodeId[0] = aNodeId2;
        aNewTriangle2.nodeId[1] = aNodeId3;
        aNewTriangle2.nodeId[2] = aNodeId4;
        aNewTriangle2.neighborId[0] = anAdvTriangle.neighborId[1];
        aNewTriangle2.neighborId[1] = aNewTriangleId3;
        aNewTriangle2.neighborId[2] = aNewTriangleId1;
        aNewTriangle2.nodeNotId[0] = anAdvTriangle.nodeNotId[1];
        aNewTriangle2.nodeNotId[1] = 1;
        aNewTriangle2.nodeNotId[2] = 0;
        aNewTriangle2.tetrahedronId = aNewElementId;
        aNewTriangle2.nodeNotId4 = aNodeId1;

        aNewTriangle2.build(m_volumeMesh);

        // Correct connectivity
        aTriangle2.neighborId[(anAdvTriangle.nodeNotId[1] + 1) % 3] = aNewTriangleId2;
        aTriangle2.nodeNotId[(anAdvTriangle.nodeNotId[1] + 1) % 3] = 2;

        // Add triangle to rtree
        this->addTriangleToRtree(aNewTriangle2, aTolerance);

        // Add triangle 3
        SAdvancingFrontTriangle aNewTriangle3;
        aNewTriangle3.id = aNewTriangleId3;
        aNewTriangle3.remove = false;
        aNewTriangle3.boundary = false;
        aNewTriangle3.nodeId[0] = aNodeId3;
        aNewTriangle3.nodeId[1] = aNodeId1;
        aNewTriangle3.nodeId[2] = aNodeId4;
        aNewTriangle3.neighborId[0] = anAdvTriangle.neighborId[2];
        aNewTriangle3.neighborId[1] = aNewTriangleId1;
        aNewTriangle3.neighborId[2] = aNewTriangleId2;
        aNewTriangle3.nodeNotId[0] = anAdvTriangle.nodeNotId[2];
        aNewTriangle3.nodeNotId[1] = 1;
        aNewTriangle3.nodeNotId[2] = 0;
        aNewTriangle3.tetrahedronId = aNewElementId;
        aNewTriangle3.nodeNotId4 = aNodeId2;

        aNewTriangle3.build(m_volumeMesh);

        // Correct connectivity
        aTriangle3.neighborId[(anAdvTriangle.nodeNotId[2] + 1) % 3] = aNewTriangleId3;
        aTriangle3.nodeNotId[(anAdvTriangle.nodeNotId[2] + 1) % 3] = 2;

        // Add triangle to rtree
        this->addTriangleToRtree(aNewTriangle3, aTolerance);

        // Remove this triangle
        this->removeTriangle(anAdvTriangle, aTolerance);

        m_anAdvFront.push_back(aNewTriangle1);
        m_anAdvFront.push_back(aNewTriangle2);
        m_anAdvFront.push_back(aNewTriangle3);

        this->cleanDuplicateTriangles(sTriangles, aTolerance);
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::generate(CMshMesh<Real>& aSurfaceMesh, const Integer maxNbElements, Real meshSize, Real minQuality, const Real aTolerance)
    {
        std::vector<CGeoCoordinate<double>> sInteriorPoints;

        ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
        aAnaFunction.set(meshSize);

        return this->generate(aSurfaceMesh, maxNbElements, sInteriorPoints, aAnaFunction, minQuality, aTolerance);
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::generate(CMshMesh<Real>& aSurfaceMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, Real meshSize, Real minQuality, const Real aTolerance)
    {
        ENigMA::analytical::CAnaFunction<Real> aAnaFunction;
        aAnaFunction.set(meshSize);

        return this->generate(aSurfaceMesh, maxNbElements, sInteriorPoints, aAnaFunction, minQuality, aTolerance);
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::generate(CMshMesh<Real>& aSurfaceMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minQuality, const Real aTolerance)
    {
        m_bStop = false;

        m_previousNbElements = 0;

        // Add boundary nodes to volume mesh
        m_volumeMesh.reset();
        m_boundingBox.reset();

        for (Integer i = 0; i < aSurfaceMesh.nbNodes(); ++i) {

            Integer aNodeId = aSurfaceMesh.nodeId(i);
            CMshNode<Real>& aNode = aSurfaceMesh.node(aNodeId);

            m_volumeMesh.addNode(aNodeId, aNode);

            // Add to bounding box
            m_boundingBox.addCoordinate(aNode);
        }

        m_boundingBox.grow(aTolerance);

        // Add boundary to advancing front and rtree
        m_anAdvFront.clear();

        m_anAdvFront.reserve(aSurfaceMesh.nbElements() * 50);

        m_tree.reset();

        std::map<Integer, Integer> newTriangleIds;

        m_nextTriangleId = 0;

        for (Integer i = 0; i < aSurfaceMesh.nbElements(); ++i) {

            Integer anElementId = aSurfaceMesh.elementId(i);

            CMshElement<Real>& anElement = aSurfaceMesh.element(anElementId);

            if (anElement.elementType() == ET_TRIANGLE)
                newTriangleIds[anElementId] = m_nextTriangleId++;
        }

        m_nextTriangleId = 0;

        for (Integer i = 0; i < aSurfaceMesh.nbElements(); ++i) {

            Integer anElementId = aSurfaceMesh.elementId(i);

            CMshElement<Real>& anElement = aSurfaceMesh.element(anElementId);

            if (anElement.elementType() == ET_TRIANGLE) {

                SAdvancingFrontTriangle anAdvTriangle;

                anAdvTriangle.id = m_nextTriangleId++;
                anAdvTriangle.remove = false;
                anAdvTriangle.boundary = true;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                    anAdvTriangle.nodeId[j] = anElement.nodeId(j);

                anAdvTriangle.tetrahedronId = std::numeric_limits<Integer>::max();
                anAdvTriangle.nodeNotId4 = std::numeric_limits<Integer>::max();

                anAdvTriangle.build(aSurfaceMesh);

                for (Integer j = 0; j < anElement.nbFaceIds(); ++j) {

                    anAdvTriangle.neighborId[j] = std::numeric_limits<Integer>::max();

                    Integer aFaceId = anElement.faceId(j);
                    CMshFace<Real>& aFace = aSurfaceMesh.face(aFaceId);

                    if (aFace.hasPair()) {
                        Integer aPairFaceId = aFace.pairFaceId();
                        Integer aNeighborId = aSurfaceMesh.face(aPairFaceId).elementId();

                        if (newTriangleIds.find(aNeighborId) != newTriangleIds.end()) {

                            anAdvTriangle.neighborId[j] = newTriangleIds.at(aNeighborId);

                            CMshElement<Real>& aNeighbor = aSurfaceMesh.element(aNeighborId);

                            for (Integer k = 0; k < aNeighbor.nbNodeIds(); ++k) {
                                if (aNeighbor.nodeId(k) != aFace.nodeId(0) && aNeighbor.nodeId(k) != aFace.nodeId(1)) {
                                    anAdvTriangle.nodeNotId[j] = k;
                                    break;
                                }
                            }

                        } else
                            std::cout << "Error: element id = " << aNeighborId << " not found!" << std::endl;

                    } else {

                        throw std::runtime_error("Boundary is open!");
                    }
                }

                m_anAdvFront.push_back(anAdvTriangle);

                this->addTriangleToRtree(anAdvTriangle, aTolerance);
            }
        }

        // Add interior nodes
        for (Integer i = 0; i < static_cast<Integer>(sInteriorPoints.size()); ++i) {

            Integer aNewNodeId = m_volumeMesh.nextNodeId();
            CMshNode<Real> aNewNode = sInteriorPoints[i];

            m_volumeMesh.addNode(aNewNodeId, aNewNode);

            SNode anInteriorNode;

            anInteriorNode.id = static_cast<Integer>(m_innerNodes.size());
            anInteriorNode.remove = false;

            anInteriorNode.nodeId = aNewNodeId;

            m_innerNodes.push_back(anInteriorNode);
        }

        // Start meshing interior

        Integer maxElem = maxNbElements;

        bool res = true;

        res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 1.00, 0.10, false, 0, aTolerance) : res = false;
        res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00, 0.50, 1.35, 0.01, false, 0, aTolerance) : res = false;

        if (res) {

            for (Integer i = 0; i < 30; ++i) {

                if (i < 10)
                    res ? res = this->repair(meshSizeFunc, 0.90 + i * 0.05, 0.0, aTolerance) : res = false;
                else
                    res ? res = this->repair(meshSizeFunc, 1.20 - i * 0.05, 0.0, aTolerance) : res = false;

                if (i < 4) {
                    res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 1.00, 0.15, false, this->getFirstIndex(), aTolerance) : res = false;
                    res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00, 1.00, 1.00, 0.15, false, 0, aTolerance) : res = false;
                } else if (i < 8) {
                    res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00 - 0.05 * i, 1.00, 1.00, 0.15, false, this->getFirstIndex(), aTolerance) : res = false;
                    res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00 - 0.05 * i, 1.00, 1.00, 0.15, false, 0, aTolerance) : res = false;
                }

                res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00, 0.50, 1.50, 0.00, false, this->getFirstIndex(), aTolerance) : res = false;
                res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00, 0.50, 1.50, 0.00, false, 0, aTolerance) : res = false;

                res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00, 0.25, 2.50, 0.00, false, this->getFirstIndex(), aTolerance) : res = false;
                res ? res = this->advancingFrontMeshing(meshSizeFunc, maxElem, 1.00, 0.25, 2.50, 0.00, false, 0, aTolerance) : res = false;

                if (!res)
                    break;

                if (this->frontSize() == 0)
                    break;
            }
        }

        m_volumeMesh.removeDanglingNodes();
        m_volumeMesh.renumber();

        if (this->frontSize() == 0) {
            m_volumeMesh.generateFaces(aTolerance);
            return true;
        } else {
            std::cout << "Meshing error!" << std::endl;
            return false;
        }
    }

    template <typename Real>
    Integer CMshTetrahedronMesher<Real>::frontSize()
    {

        Integer n = 0;

        for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {
            if (!m_anAdvFront[i].remove)
                n++;
        }

        return n;
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::advancingFrontMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real sizeFactor, Real shrinkFactor, Real expandFactor, Real minQuality, const bool bCheckDelaunay, Integer firstIndex, const Real aTolerance)
    {
        bool res = true;

        Real x, y, z;

        meshSizeFunc.removeAllVariables();

        meshSizeFunc.defineVariable("x", x);
        meshSizeFunc.defineVariable("y", y);
        meshSizeFunc.defineVariable("z", z);

        Real sumMeshSize = 0;
        Real averageMeshSize = 0;
        Integer nMeshSize = 0;

        for (Integer i = firstIndex; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {

            this->checkUpdate();

            if (m_bStop)
                return false;

            if (maxNbElements > 0) {
                if (m_volumeMesh.nbElements() >= maxNbElements) {
                    std::cout << "Max number of elements (" << maxNbElements << ") reached!" << std::endl;
                    return false;
                }
            }

            Integer anAdvTriangleId = i;

            if (m_anAdvFront[anAdvTriangleId].remove)
                continue;

            SAdvancingFrontTriangle& anAdvTriangle = m_anAdvFront[anAdvTriangleId];

            if (anAdvTriangle.neighborId[0] == std::numeric_limits<Integer>::max() || anAdvTriangle.neighborId[1] == std::numeric_limits<Integer>::max() || anAdvTriangle.neighborId[2] == std::numeric_limits<Integer>::max())
                continue;

            Integer aNodeId1 = anAdvTriangle.nodeId[0];
            Integer aNodeId2 = anAdvTriangle.nodeId[1];
            Integer aNodeId3 = anAdvTriangle.nodeId[2];

            CMshNode<Real>& aNode1 = m_volumeMesh.node(aNodeId1);
            CMshNode<Real>& aNode2 = m_volumeMesh.node(aNodeId2);
            CMshNode<Real>& aNode3 = m_volumeMesh.node(aNodeId3);

            Integer aNodeId4 = m_anAdvFront[anAdvTriangle.neighborId[0]].nodeId[anAdvTriangle.nodeNotId[0]];
            Integer aNodeId5 = m_anAdvFront[anAdvTriangle.neighborId[1]].nodeId[anAdvTriangle.nodeNotId[1]];
            Integer aNodeId6 = m_anAdvFront[anAdvTriangle.neighborId[2]].nodeId[anAdvTriangle.nodeNotId[2]];

            CMshNode<Real>& aNode4 = m_volumeMesh.node(aNodeId4);

            CMshNode<Real> aMidNode = (aNode1 + aNode2 + aNode3) / 3.0;

            x = aMidNode.x();
            y = aMidNode.y();
            z = aMidNode.z();

            // Use inverse of normal vector
            CGeoVector<Real> v = -((aNode2 - aNode1).cross(aNode3 - aNode1));

            Real localMeshSize = std::max(meshSizeFunc.evaluate(), static_cast<Real>(v.norm() * 0.7));

            sumMeshSize += localMeshSize;
            nMeshSize++;
            averageMeshSize = sumMeshSize / nMeshSize;

            if (averageMeshSize > localMeshSize)
                localMeshSize = std::min(averageMeshSize, static_cast<Real>(2.0 * localMeshSize));
            else
                localMeshSize = std::max(averageMeshSize, static_cast<Real>(0.5 * localMeshSize));

            v.normalize();

            // Regular tetrahedron (height to edge ratio)
            v *= sqrt(2.0 / 3.0);

            // Add point to form tetrahedra with correct spacing
            CMshNode<Real> aNewNode = aMidNode + v * localMeshSize * sizeFactor;
            Integer aNewNodeId = m_volumeMesh.nextNodeId();

            // Get closest triangles
            CMshNode<Real> anAuxNode = aMidNode + v * localMeshSize * sizeFactor * 1.5;

            CGeoBoundingBox<Real> aBoundingBox;
            aBoundingBox.addCoordinate(aNode1);
            aBoundingBox.addCoordinate(aNode2);
            aBoundingBox.addCoordinate(aNode3);
            aBoundingBox.addCoordinate(anAuxNode);
            aBoundingBox.grow(aTolerance);

            std::vector<Integer> sTriangles;
            m_tree.find(sTriangles, aBoundingBox);

            sTriangles.erase(std::remove(sTriangles.begin(), sTriangles.end(), anAdvTriangleId), sTriangles.end());

            // Check if a node exists in proximity
            std::vector<Integer> sNodes;
            this->findClosestNodes(sTriangles, sNodes);

            // Meshing priority
            // Priority = 1: close hole
            // Priority = 2: other nodes in vicinity (4, 5, 6, other node)
            // Priority = 3: new node forming correct spacing

            if (aNodeId4 == aNodeId5 && aNodeId5 == aNodeId6) {

                CMshTetrahedron<Real> aNewTetrahedron;

                aNewTetrahedron.addVertex(aNode1);
                aNewTetrahedron.addVertex(aNode2);
                aNewTetrahedron.addVertex(aNode3);
                aNewTetrahedron.addVertex(aNode4);

                aNewTetrahedron.calculateVolume();

                if (aNewTetrahedron.volume() > aTolerance * aTolerance * aTolerance) {

                    Integer wNodeId;

                    if (this->triangleOk(anAdvTriangle, aNode1, aNode2, aNode4, sTriangles, aTolerance) && this->triangleOk(anAdvTriangle, aNode2, aNode3, aNode4, sTriangles, aTolerance) && this->triangleOk(anAdvTriangle, aNode3, aNode1, aNode4, sTriangles, aTolerance) && !this->tetrahedronContainsNode(aNode1, aNode2, aNode3, aNode4, wNodeId, sNodes, aTolerance)) {

                        this->addTetrahedron(anAdvTriangle, aNodeId4, sTriangles, aTolerance);
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

                if (aNodeId == aNodeId1 || aNodeId == aNodeId2 || aNodeId == aNodeId3)
                    continue;

                CMshNode<Real>& aNode = m_volumeMesh.node(aNodeId);

                Real d = (aNode - aNewNode).norm();

                // Use closest node
                if (d < localMeshSize * sizeFactor * expandFactor) {

                    CMshTetrahedron<Real> aNewTetrahedron1;

                    aNewTetrahedron1.addVertex(aNode1);
                    aNewTetrahedron1.addVertex(aNode2);
                    aNewTetrahedron1.addVertex(aNode3);
                    aNewTetrahedron1.addVertex(aNode);

                    aNewTetrahedron1.calculateVolume();

                    if (aNewTetrahedron1.volume() > aTolerance * aTolerance * aTolerance) {

                        aNewTetrahedron1.calculateQuality();

                        Real q1 = aNewTetrahedron1.quality();

                        if (q1 > qmax)
                            q1 += this->triangleOk(anAdvTriangle, aNode1, aNode2, aNode, sTriangles, aTolerance) ? 0.0 : -2.0;

                        if (q1 > qmax)
                            q1 += this->triangleOk(anAdvTriangle, aNode2, aNode3, aNode, sTriangles, aTolerance) ? 0.0 : -2.0;

                        if (q1 > qmax)
                            q1 += this->triangleOk(anAdvTriangle, aNode3, aNode1, aNode, sTriangles, aTolerance) ? 0.0 : -2.0;

                        if (q1 > qmax) {
                            Integer wNodeId;
                            q1 += !this->tetrahedronContainsNode(aNode1, aNode2, aNode3, aNode, wNodeId, sNodes, aTolerance) ? 0.0 : -2.0;
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

                for (Integer k = 0; k < static_cast<Integer>(m_innerNodes.size()); ++k) {

                    if (m_innerNodes[k].remove)
                        continue;

                    if (m_innerNodes[k].nodeId == anExistingNodeId)
                        m_innerNodes[k].remove = true;
                }

                this->addTetrahedron(anAdvTriangle, anExistingNodeId, sTriangles, aTolerance);
                res = true;

            } else {

                CGeoLine<Real> aLine(aMidNode, aNewNode);

                Real dmin = findShortestDistance(sTriangles, aLine, anAdvTriangleId, aTolerance);

                if (dmin > localMeshSize * sizeFactor * shrinkFactor * 0.25 && dmin < localMeshSize * sizeFactor * expandFactor)
                    aNewNode = aMidNode + v * dmin * 0.5;

                CMshTetrahedron<Real> aNewTetrahedron2;

                aNewTetrahedron2.addVertex(aNode1);
                aNewTetrahedron2.addVertex(aNode2);
                aNewTetrahedron2.addVertex(aNode3);
                aNewTetrahedron2.addVertex(aNewNode);

                aNewTetrahedron2.calculateVolume();

                if (aNewTetrahedron2.volume() > aTolerance * aTolerance * aTolerance) {

                    aNewTetrahedron2.calculateQuality();

                    Real q2 = aNewTetrahedron2.quality();

                    if (q2 > 1.5 * minQuality)
                        q2 += this->triangleOk(anAdvTriangle, aNode1, aNode2, aNewNode, sTriangles, aTolerance) ? 0.0 : -2.0;

                    if (q2 > 1.5 * minQuality)
                        q2 += this->triangleOk(anAdvTriangle, aNode2, aNode3, aNewNode, sTriangles, aTolerance) ? 0.0 : -2.0;

                    if (q2 > 1.5 * minQuality)
                        q2 += this->triangleOk(anAdvTriangle, aNode3, aNode1, aNewNode, sTriangles, aTolerance) ? 0.0 : -2.0;

                    if (q2 > 1.5 * minQuality) {
                        Integer wNodeId;
                        q2 += !this->tetrahedronContainsNode(aNode1, aNode2, aNode3, aNewNode, wNodeId, sNodes, aTolerance) ? 0.0 : -2.0;
                    }

                    if (q2 > 1.5 * minQuality && bCheckDelaunay)
                        q2 += this->checkDelaunay(aNewNode, aTolerance) ? 0.0 : -2.0;

                    if (q2 > 1.5 * minQuality) {

                        // Create a new node
                        m_volumeMesh.addNode(aNewNodeId, aNewNode);

                        this->addTetrahedron(anAdvTriangle, aNewNodeId, sTriangles, aTolerance);
                        res = true;

                        if (!m_boundingBox.contains(aNewNode, aTolerance)) {
                            throw std::runtime_error("Node is outside boundary!");
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
    bool CMshTetrahedronMesher<Real>::rebuildConnectivity(const Real aTolerance)
    {
        CGeoBoundingBox<Real> aBoundingBox;

        // Build connectivity
        for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {

            Integer anAdvTriangleId = i;

            SAdvancingFrontTriangle& anAdvTriangle = m_anAdvFront[anAdvTriangleId];

            if (anAdvTriangle.remove)
                continue;

            Integer aNodeId1 = anAdvTriangle.nodeId[0];
            Integer aNodeId2 = anAdvTriangle.nodeId[1];
            Integer aNodeId3 = anAdvTriangle.nodeId[2];

            CMshNode<Real>& aNode1 = m_volumeMesh.node(aNodeId1);
            CMshNode<Real>& aNode2 = m_volumeMesh.node(aNodeId2);
            CMshNode<Real>& aNode3 = m_volumeMesh.node(aNodeId3);

            aBoundingBox.reset();

            aBoundingBox.addCoordinate(aNode1);
            aBoundingBox.addCoordinate(aNode2);
            aBoundingBox.addCoordinate(aNode3);

            aBoundingBox.grow(aTolerance);

            std::vector<Integer> sTriangles;
            m_tree.find(sTriangles, aBoundingBox);

            for (Integer j = 0; j < static_cast<Integer>(sTriangles.size()); ++j) {

                Integer anotherTriangleId = sTriangles[j];

                // Discard same triangle
                if (anotherTriangleId == anAdvTriangle.id)
                    continue;

                SAdvancingFrontTriangle& anotherTriangle = m_anAdvFront[anotherTriangleId];

                if (anotherTriangle.remove)
                    continue;

                this->pairEdges(anAdvTriangle, anotherTriangle);
            }

            if (anAdvTriangle.neighborId[0] == std::numeric_limits<Integer>::max() || anAdvTriangle.neighborId[1] == std::numeric_limits<Integer>::max() || anAdvTriangle.neighborId[2] == std::numeric_limits<Integer>::max()) {
                throw std::runtime_error("Boundary is open!");
            }
        }
        
        return true;
    }

    template <typename Real>
    Integer CMshTetrahedronMesher<Real>::getFirstIndex()
    {

        Integer aFirstIndex = 0;

        Real minArea = std::numeric_limits<Real>::max();

        for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {

            SAdvancingFrontTriangle& anAdvTriangle = m_anAdvFront[i];

            if (anAdvTriangle.remove)
                continue;

            CGeoTriangle<Real>& aTriangle = anAdvTriangle.triangle;

            aTriangle.calculateArea();

            if (aTriangle.area() < minArea) {
                minArea = aTriangle.area();
                aFirstIndex = i;
            }
        }

        return aFirstIndex;
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::reduceFront(const Real aTolerance)
    {

        std::vector<SAdvancingFrontTriangle> aReducedAdvFront;

        std::map<Integer, Integer> aAdvFrontMapId;

        Integer aNewAdvFrontId = 0;

        for (Integer i = 0; i < m_anAdvFront.size(); ++i) {

            SAdvancingFrontTriangle& anAdvTriangle = m_anAdvFront[i];

            if (anAdvTriangle.remove)
                continue;

            aAdvFrontMapId[i] = aNewAdvFrontId;
            aNewAdvFrontId++;

            aReducedAdvFront.push_back(m_anAdvFront[i]);
        }

        m_tree.reset();

        m_nextTriangleId = aReducedAdvFront.size();

        m_anAdvFront.clear();

        for (Integer i = 0; i < aReducedAdvFront.size(); ++i) {

            SAdvancingFrontTriangle& anAdvTriangle = aReducedAdvFront[i];

            anAdvTriangle.id = aAdvFrontMapId.at(anAdvTriangle.id);
            anAdvTriangle.neighborId[0] = aAdvFrontMapId.at(anAdvTriangle.neighborId[0]);
            anAdvTriangle.neighborId[1] = aAdvFrontMapId.at(anAdvTriangle.neighborId[1]);
            anAdvTriangle.neighborId[2] = aAdvFrontMapId.at(anAdvTriangle.neighborId[2]);

            m_anAdvFront.push_back(anAdvTriangle);

            this->addTriangleToRtree(anAdvTriangle, aTolerance);
        }
    }

    template <typename Real>
    bool CMshTetrahedronMesher<Real>::repair(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, const Real sizeFactor, const Real minQuality, const Real aTolerance)
    {

        Real x, y, z;

        meshSizeFunc.removeAllVariables();

        meshSizeFunc.defineVariable("x", x);
        meshSizeFunc.defineVariable("y", y);
        meshSizeFunc.defineVariable("z", z);

        std::vector<CGeoCoordinate<Real>> sSteinerPoints;

        std::vector<Integer> sElementsToRemove;

        for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {

            this->checkUpdate();

            if (m_bStop)
                return false;

            SAdvancingFrontTriangle& anAdvTriangle = m_anAdvFront[i];

            if (anAdvTriangle.remove)
                continue;

            if (anAdvTriangle.tetrahedronId == std::numeric_limits<Integer>::max())
                continue;

            for (Integer j = i + 1; j < static_cast<Integer>(m_anAdvFront.size()); ++j) {

                SAdvancingFrontTriangle& anotherAdvTriangle = m_anAdvFront[j];

                if (anotherAdvTriangle.remove)
                    continue;

                if (anotherAdvTriangle.tetrahedronId == std::numeric_limits<Integer>::max())
                    continue;

                if (std::find(sElementsToRemove.begin(), sElementsToRemove.end(), anAdvTriangle.tetrahedronId) == sElementsToRemove.end())
                    sElementsToRemove.push_back(anAdvTriangle.tetrahedronId);

                if (std::find(sElementsToRemove.begin(), sElementsToRemove.end(), anotherAdvTriangle.tetrahedronId) == sElementsToRemove.end())
                    sElementsToRemove.push_back(anotherAdvTriangle.tetrahedronId);

                for (Integer k = 0; k < 3; ++k) {

                    Integer aNodeId1 = anAdvTriangle.nodeId[(k + 0) % 3];
                    Integer aNodeId2 = anAdvTriangle.nodeId[(k + 1) % 3];

                    CMshNode<Real>& aNode1 = m_volumeMesh.node(aNodeId1);
                    CMshNode<Real>& aNode2 = m_volumeMesh.node(aNodeId2);

                    CGeoLine<Real> aLine1(aNode1, aNode2);

                    for (Integer l = 0; l < 3; l++) {

                        Integer aNodeId3 = anotherAdvTriangle.nodeId[(l + 0) % 3];
                        Integer aNodeId4 = anotherAdvTriangle.nodeId[(l + 1) % 3];

                        if (aNodeId1 == aNodeId3 || aNodeId1 == aNodeId4 || aNodeId2 == aNodeId3 || aNodeId2 == aNodeId4)
                            continue;

                        CMshNode<Real>& aNode3 = m_volumeMesh.node(aNodeId3);
                        CMshNode<Real>& aNode4 = m_volumeMesh.node(aNodeId4);

                        CGeoLine<Real> aLine2(aNode3, aNode4);

                        CGeoCoordinate<Real> aPoint5;
                        CGeoCoordinate<Real> aPoint6;

                        Real aDistance = 0.0;

                        if (aLine1.distance(aLine2, aPoint5, aPoint6, aDistance, aTolerance)) {

                            aLine1.calculateLength();
                            aLine2.calculateLength();

                            Real aLength = (aLine1.length() + aLine2.length()) * 0.5;

                            if (aDistance < aLength * 0.20) {

                                CGeoCoordinate<Real> aNewPoint = (aPoint6 + aPoint5) * 0.5;

                                bool bAddNode = true;

                                for (Integer m = 0; m < static_cast<Integer>(sSteinerPoints.size()); m++) {

                                    if ((sSteinerPoints[m] - aNewPoint).norm() < aLength * 0.50) {
                                        bAddNode = false;
                                        break;
                                    }
                                }

                                if (bAddNode)
                                    sSteinerPoints.push_back(aNewPoint);
                            }
                        }
                    }
                }
            }
        }

        for (Integer i = 0; i < static_cast<Integer>(sSteinerPoints.size()); ++i) {

            this->checkUpdate();

            if (m_bStop)
                return false;

            x = sSteinerPoints[i].x();
            y = sSteinerPoints[i].y();
            z = sSteinerPoints[i].z();

            Real localMeshSize = meshSizeFunc.evaluate();

            CGeoSphere<Real> aSphere(sSteinerPoints[i], localMeshSize * sizeFactor);

            for (Integer j = 0; j < m_volumeMesh.nbElements(); ++j) {

                Integer anElementId = m_volumeMesh.elementId(j);

                CMshElement<Real>& anElement = m_volumeMesh.element(anElementId);

                for (Integer k = 0; k < anElement.nbNodeIds(); ++k) {

                    Integer aNodeId = anElement.nodeId(k);

                    CMshNode<Real>& aNode = m_volumeMesh.node(aNodeId);

                    if (aSphere.contains(aNode, aTolerance)) {

                        if (std::find(sElementsToRemove.begin(), sElementsToRemove.end(), anElementId) == sElementsToRemove.end())
                            sElementsToRemove.push_back(anElementId);
                    }
                }
            }
        }

        for (Integer i = 0; i < static_cast<Integer>(sElementsToRemove.size()); ++i) {

            this->checkUpdate();

            if (m_bStop)
                return false;

            Integer anElementId = sElementsToRemove[i];

            std::vector<CMshFace<Real>> sFaces;

            CMshElement<Real>& anElement = m_volumeMesh.element(anElementId);

            anElement.generateFaces(sFaces);

            for (Integer j = 0; j < static_cast<Integer>(sFaces.size()); ++j) {

                SAdvancingFrontTriangle anAdvTriangle;

                anAdvTriangle.id = m_nextTriangleId++;
                anAdvTriangle.remove = false;
                anAdvTriangle.boundary = false;

                anAdvTriangle.nodeId[0] = sFaces[j].nodeId(0);
                anAdvTriangle.nodeId[1] = sFaces[j].nodeId(1);
                anAdvTriangle.nodeId[2] = sFaces[j].nodeId(2);

                anAdvTriangle.neighborId[0] = std::numeric_limits<Integer>::max();
                anAdvTriangle.neighborId[1] = std::numeric_limits<Integer>::max();
                anAdvTriangle.neighborId[2] = std::numeric_limits<Integer>::max();

                anAdvTriangle.tetrahedronId = std::numeric_limits<Integer>::max();
                anAdvTriangle.nodeNotId4 = std::numeric_limits<Integer>::max();

                anAdvTriangle.build(m_volumeMesh);

                m_anAdvFront.push_back(anAdvTriangle);

                this->addTriangleToRtree(anAdvTriangle, aTolerance);
            }
        }

        for (Integer i = 0; i < static_cast<Integer>(m_anAdvFront.size()); ++i) {

            this->checkUpdate();

            if (m_bStop)
                return false;

            SAdvancingFrontTriangle& anAdvTriangle = m_anAdvFront[i];

            if (anAdvTriangle.remove)
                continue;

            std::vector<Integer> sNodeIds;

            sNodeIds.push_back(anAdvTriangle.nodeId[0]);
            sNodeIds.push_back(anAdvTriangle.nodeId[1]);
            sNodeIds.push_back(anAdvTriangle.nodeId[2]);

            std::sort(sNodeIds.begin(), sNodeIds.end());

            for (Integer j = i + 1; j < static_cast<Integer>(m_anAdvFront.size()); ++j) {

                SAdvancingFrontTriangle& anotherAdvTriangle = m_anAdvFront[j];

                if (anotherAdvTriangle.remove)
                    continue;

                std::vector<Integer> sOtherNodeIds;

                sOtherNodeIds.push_back(anotherAdvTriangle.nodeId[0]);
                sOtherNodeIds.push_back(anotherAdvTriangle.nodeId[1]);
                sOtherNodeIds.push_back(anotherAdvTriangle.nodeId[2]);

                std::sort(sOtherNodeIds.begin(), sOtherNodeIds.end());

                if (sNodeIds == sOtherNodeIds) {

                    this->removeTriangle(anAdvTriangle, aTolerance);
                    this->removeTriangle(anotherAdvTriangle, aTolerance);
                }
            }
        }

        this->rebuildConnectivity(aTolerance);

        for (Integer i = 0; i < static_cast<Integer>(sElementsToRemove.size()); ++i)
            m_volumeMesh.removeElement(sElementsToRemove[i]);

        return true;
    }

    template <typename Real>
    CMshMesh<Real>& CMshTetrahedronMesher<Real>::mesh()
    {

        return m_volumeMesh;
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::setIntervals(const Integer timeInterval, const Integer dataInterval)
    {

        m_begin = clock();

        m_timeInterval = timeInterval;
        m_dataInterval = dataInterval;
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::stopMeshing()
    {

        m_bStop = true;
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::relaxNodes(const Real aFactor, const Real aTolerance)
    {

        std::map<Integer, bool> bBoundaryNode;

        for (Integer i = 0; i < m_volumeMesh.nbNodes(); ++i) {

            Integer aNodeId = m_volumeMesh.nodeId(i);
            bBoundaryNode[aNodeId] = false;
        }

        for (Integer i = 0; i < m_volumeMesh.nbFaces(); ++i) {

            Integer aFaceId = m_volumeMesh.faceId(i);
            CMshFace<Real>& aFace = m_volumeMesh.face(aFaceId);

            if (aFace.hasPair())
                continue;

            for (Integer j = 0; j < aFace.nbNodeIds(); ++j) {

                Integer aNodeId = aFace.nodeId(j);
                bBoundaryNode[aNodeId] = true;
            }
        }

        CMshMeshQuery<Real> aMeshQuery(m_volumeMesh);

        std::vector<Integer> sElementIds;

        for (Integer i = 0; i < m_volumeMesh.nbNodes(); ++i) {

            Integer aMovingNodeId = m_volumeMesh.nodeId(i);
            CMshNode<Real>& aMovingNode = m_volumeMesh.node(aMovingNodeId);

            if (bBoundaryNode.at(aMovingNodeId))
                continue;

            aMeshQuery.elementsSharingNode(aMovingNodeId, sElementIds);

            std::vector<Integer> sNodeIds;
            sNodeIds.push_back(aMovingNodeId);

            Real aTotalVolume = 0.0;

            CMshNode<Real> aNewNode(0.0, 0.0, 0.0);

            for (Integer j = 0; j < static_cast<Integer>(sElementIds.size()); ++j) {

                Integer anElementId = sElementIds[j];

                CMshElement<Real>& anElement = m_volumeMesh.element(anElementId);

                if (anElement.elementType() != ET_TETRAHEDRON)
                    continue;

                CMshNode<Real> aNode1 = m_volumeMesh.node(anElement.nodeId(0));
                CMshNode<Real> aNode2 = m_volumeMesh.node(anElement.nodeId(1));
                CMshNode<Real> aNode3 = m_volumeMesh.node(anElement.nodeId(2));
                CMshNode<Real> aNode4 = m_volumeMesh.node(anElement.nodeId(3));

                CMshTetrahedron<Real> aTetrahedron;

                aTetrahedron.addVertex(aNode1);
                aTetrahedron.addVertex(aNode2);
                aTetrahedron.addVertex(aNode3);
                aTetrahedron.addVertex(aNode4);

                aTetrahedron.calculateVolume();

                aNewNode += aNode1;
                aNewNode += aNode2;
                aNewNode += aNode3;
                aNewNode += aNode4;
                aNewNode /= 4.0;

                aNewNode *= aTetrahedron.volume();

                aTotalVolume += aTetrahedron.volume();
            }

            if (aTotalVolume > 0.0)
                aNewNode /= aTotalVolume;

            CGeoVector<Real> v = (aNewNode - aMovingNode);

            std::vector<CMshNode<Real>> sNewNodes;
            sNewNodes.push_back(aMovingNode + v * aFactor);
            sNewNodes.push_back(aMovingNode - v * aFactor);

            for (Integer n = 0; n < static_cast<Integer>(sNewNodes.size()); n++) {

                aNewNode = sNewNodes[n];

                Real minq1 = 1.0;
                Real minq2 = 1.0;

                bool bMoveNode = true;

                for (Integer k = 0; k < static_cast<Integer>(sElementIds.size()); ++k) {

                    CMshElement<Real>& aModElement = m_volumeMesh.element(sElementIds[k]);

                    if (aModElement.elementType() != ET_TETRAHEDRON)
                        continue;

                    Integer aNodeId1 = aModElement.nodeId(0);
                    Integer aNodeId2 = aModElement.nodeId(1);
                    Integer aNodeId3 = aModElement.nodeId(2);
                    Integer aNodeId4 = aModElement.nodeId(3);

                    CMshNode<Real>& aNode1 = m_volumeMesh.node(aNodeId1);
                    CMshNode<Real>& aNode2 = m_volumeMesh.node(aNodeId2);
                    CMshNode<Real>& aNode3 = m_volumeMesh.node(aNodeId3);
                    CMshNode<Real>& aNode4 = m_volumeMesh.node(aNodeId4);

                    CGeoCoordinate<Real> aVertex1 = aNode1;
                    CGeoCoordinate<Real> aVertex2 = aNode2;
                    CGeoCoordinate<Real> aVertex3 = aNode3;
                    CGeoCoordinate<Real> aVertex4 = aNode4;

                    CMshTetrahedron<Real> aTetrahedron1;

                    aTetrahedron1.addVertex(aVertex1);
                    aTetrahedron1.addVertex(aVertex2);
                    aTetrahedron1.addVertex(aVertex3);
                    aTetrahedron1.addVertex(aVertex4);

                    aTetrahedron1.calculateQuality();

                    Real q1 = aTetrahedron1.quality();

                    if (aNodeId1 == aMovingNodeId)
                        aVertex1 = aNewNode;

                    if (aNodeId2 == aMovingNodeId)
                        aVertex2 = aNewNode;

                    if (aNodeId3 == aMovingNodeId)
                        aVertex3 = aNewNode;

                    if (aNodeId4 == aMovingNodeId)
                        aVertex4 = aNewNode;

                    CMshTetrahedron<Real> aTetrahedron2;

                    aTetrahedron2.addVertex(aVertex1);
                    aTetrahedron2.addVertex(aVertex2);
                    aTetrahedron2.addVertex(aVertex3);
                    aTetrahedron2.addVertex(aVertex4);

                    aTetrahedron2.calculateVolume();
                    aTetrahedron2.calculateQuality();

                    Real q2 = aTetrahedron2.quality();

                    minq1 = std::min(minq1, q1);
                    minq2 = std::min(minq2, q2);

                    if (aTetrahedron2.volume() < aTolerance * aTolerance * aTolerance) {
                        bMoveNode = false;
                        break;
                    }
                }

                if (bMoveNode && minq2 >= minq1) {
                    aMovingNode = aNewNode;
                    break;
                }
            }
        }
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::flipEdges23(const Real aTolerance)
    {

        std::map<Integer, bool> sFlipped;

        for (Integer i = 0; i < m_volumeMesh.nbElements(); ++i) {

            Integer anElementId = m_volumeMesh.elementId(i);
            sFlipped[anElementId] = false;
        }

        for (Integer i = 0; i < m_volumeMesh.nbFaces(); ++i) {

            Integer aFaceId = m_volumeMesh.faceId(i);
            CMshFace<Real> aFace = m_volumeMesh.face(aFaceId);

            if (aFace.faceType() != FT_TRIANGLE)
                continue;

            Integer anElementId = aFace.elementId();

            if (sFlipped[anElementId])
                continue;

            CMshElement<Real>& anElement = m_volumeMesh.element(anElementId);

            if (anElement.elementType() != ET_TETRAHEDRON)
                continue;

            Integer aNodeId1 = aFace.nodeId(0);
            Integer aNodeId2 = aFace.nodeId(1);
            Integer aNodeId3 = aFace.nodeId(2);

            Integer aNodeId4;

            bool bNodeId4 = false;

            for (Integer j = 0; j < anElement.nbNodeIds(); ++j) {

                if (anElement.nodeId(j) != aNodeId1 && anElement.nodeId(j) != aNodeId2 && anElement.nodeId(j) != aNodeId3) {
                    aNodeId4 = anElement.nodeId(j);
                    bNodeId4 = true;
                    break;
                }
            }

            if (!bNodeId4)
                continue;

            if (aFace.hasPair()) {

                Integer pairFaceId = aFace.pairFaceId();
                CMshFace<Real> aPairFace = m_volumeMesh.face(pairFaceId);

                Integer aNeighborId = aPairFace.elementId();

                if (sFlipped[aNeighborId])
                    continue;

                CMshElement<Real>& aNeighbor = m_volumeMesh.element(aNeighborId);

                if (aNeighbor.elementType() != ET_TETRAHEDRON)
                    continue;

                Integer aNodeId5;

                bool bNodeId5 = false;

                for (Integer j = 0; j < aNeighbor.nbNodeIds(); ++j) {

                    if (aNeighbor.nodeId(j) != aNodeId1 && aNeighbor.nodeId(j) != aNodeId2 && aNeighbor.nodeId(j) != aNodeId3) {
                        aNodeId5 = aNeighbor.nodeId(j);
                        bNodeId5 = true;
                        break;
                    }
                }

                if (!bNodeId5)
                    continue;

                CMshNode<Real>& aNode1 = m_volumeMesh.node(aNodeId1);
                CMshNode<Real>& aNode2 = m_volumeMesh.node(aNodeId2);
                CMshNode<Real>& aNode3 = m_volumeMesh.node(aNodeId3);
                CMshNode<Real>& aNode4 = m_volumeMesh.node(aNodeId4);
                CMshNode<Real>& aNode5 = m_volumeMesh.node(aNodeId5);

                // Original
                CMshTetrahedron<Real> aTetrahedron1;
                aTetrahedron1.addVertex(aNode1);
                aTetrahedron1.addVertex(aNode2);
                aTetrahedron1.addVertex(aNode3);
                aTetrahedron1.addVertex(aNode4);

                aTetrahedron1.calculateQuality();
                Real q1 = aTetrahedron1.quality();

                CMshTetrahedron<Real> aTetrahedron2;
                aTetrahedron2.addVertex(aNode3);
                aTetrahedron2.addVertex(aNode2);
                aTetrahedron2.addVertex(aNode1);
                aTetrahedron2.addVertex(aNode5);

                aTetrahedron2.calculateQuality();
                Real q2 = aTetrahedron2.quality();

                // Modified
                CMshTetrahedron<Real> aTetrahedron3;
                aTetrahedron3.addVertex(aNode1);
                aTetrahedron3.addVertex(aNode2);
                aTetrahedron3.addVertex(aNode5);
                aTetrahedron3.addVertex(aNode4);

                aTetrahedron3.calculateVolume();
                aTetrahedron3.calculateQuality();
                Real q3 = aTetrahedron3.quality();

                CMshTetrahedron<Real> aTetrahedron4;
                aTetrahedron4.addVertex(aNode2);
                aTetrahedron4.addVertex(aNode3);
                aTetrahedron4.addVertex(aNode5);
                aTetrahedron4.addVertex(aNode4);

                aTetrahedron4.calculateVolume();
                aTetrahedron4.calculateQuality();
                Real q4 = aTetrahedron4.quality();

                CMshTetrahedron<Real> aTetrahedron5;
                aTetrahedron5.addVertex(aNode1);
                aTetrahedron5.addVertex(aNode3);
                aTetrahedron5.addVertex(aNode4);
                aTetrahedron5.addVertex(aNode5);

                aTetrahedron5.calculateVolume();
                aTetrahedron5.calculateQuality();
                Real q5 = aTetrahedron5.quality();

                if (std::min(q3, std::min(q4, q5)) > std::min(q1, q2) && aTetrahedron3.volume() > aTolerance * aTolerance * aTolerance && aTetrahedron4.volume() > aTolerance * aTolerance * aTolerance && aTetrahedron5.volume() > aTolerance * aTolerance * aTolerance) {

                    // Do flip
                    m_volumeMesh.element(anElementId).setNodeId(0, aNodeId1);
                    m_volumeMesh.element(anElementId).setNodeId(1, aNodeId2);
                    m_volumeMesh.element(anElementId).setNodeId(2, aNodeId5);
                    m_volumeMesh.element(anElementId).setNodeId(3, aNodeId4);

                    m_volumeMesh.element(aNeighborId).setNodeId(0, aNodeId2);
                    m_volumeMesh.element(aNeighborId).setNodeId(1, aNodeId3);
                    m_volumeMesh.element(aNeighborId).setNodeId(2, aNodeId5);
                    m_volumeMesh.element(aNeighborId).setNodeId(3, aNodeId4);

                    // Create one more and add it
                    CMshElement<Real> aNewElement;

                    aNewElement.setElementType(ET_TETRAHEDRON);

                    aNewElement.addNodeId(aNodeId1);
                    aNewElement.addNodeId(aNodeId3);
                    aNewElement.addNodeId(aNodeId4);
                    aNewElement.addNodeId(aNodeId5);

                    Integer aNewElementId = m_volumeMesh.nextElementId();

                    m_volumeMesh.addElement(aNewElementId, aNewElement);

                    sFlipped[anElementId] = true;
                    sFlipped[aNeighborId] = true;
                    sFlipped[aNewElementId] = true;
                }
            }
        }

        m_volumeMesh.generateFaces(aTolerance);
    }

    template <typename Real>
    void CMshTetrahedronMesher<Real>::flipEdges32(const Real aTolerance)
    {

        std::map<Integer, bool> sFlipped;
        std::map<Integer, bool> bDeleteElement;

        std::map<Integer, bool> bBoundaryNode;

        for (Integer i = 0; i < m_volumeMesh.nbNodes(); ++i) {

            Integer aNodeId = m_volumeMesh.nodeId(i);
            bBoundaryNode[aNodeId] = false;
        }

        for (Integer i = 0; i < m_volumeMesh.nbFaces(); ++i) {

            Integer aFaceId = m_volumeMesh.faceId(i);
            CMshFace<Real>& aFace = m_volumeMesh.face(aFaceId);

            if (aFace.hasPair())
                continue;

            for (Integer j = 0; j < aFace.nbNodeIds(); ++j) {

                Integer aNodeId = aFace.nodeId(j);
                bBoundaryNode[aNodeId] = true;
            }
        }

        for (Integer i = 0; i < m_volumeMesh.nbElements(); ++i) {

            Integer anElementId = m_volumeMesh.elementId(i);
            sFlipped[anElementId] = false;

            bDeleteElement[anElementId] = false;
        }

        CMshMeshQuery<Real> aMeshQuery(m_volumeMesh);

        for (Integer i = 0; i < m_volumeMesh.nbElements(); ++i) {

            Integer anElementId = m_volumeMesh.elementId(i);
            CMshElement<Real>& anElement = m_volumeMesh.element(anElementId);

            if (anElement.elementType() != ET_TETRAHEDRON)
                continue;

            if (sFlipped.at(anElementId))
                continue;

            Integer aNodeId1 = anElement.nodeId(0);
            Integer aNodeId2 = anElement.nodeId(1);
            Integer aNodeId3 = anElement.nodeId(2);
            Integer aNodeId4 = anElement.nodeId(3);

            for (Integer j = 0; j < anElement.nbFaceIds(); ++j) {

                if (sFlipped.at(anElementId))
                    continue;

                Integer aFaceId = anElement.faceId(j);
                CMshFace<Real>& aFace = m_volumeMesh.face(aFaceId);

                if (aFace.faceType() != FT_TRIANGLE)
                    continue;

                std::vector<Integer> sNodeIds;
                std::vector<Integer> sElementIds;

                for (Integer k = 0; k < aFace.nbNodeIds(); ++k) {

                    if (sFlipped.at(anElementId))
                        continue;

                    if (bBoundaryNode.at(aFace.nodeId((k + 0) % 3)) && bBoundaryNode.at(aFace.nodeId((k + 1) % 3)))
                        continue;

                    aMeshQuery.elementsSharingNodes(aFace.nodeId((k + 0) % 3), aFace.nodeId((k + 1) % 3), sElementIds);

                    if (sElementIds.size() == 3) {

                        Integer aNodeId5;

                        bool bNodeFound = false;

                        for (Integer l = 0; l < static_cast<Integer>(sElementIds.size()); l++) {

                            Integer anOtherElementId = sElementIds[l];
                            CMshElement<Real>& anOtherElement = m_volumeMesh.element(anOtherElementId);

                            if (anOtherElement.elementType() != ET_TETRAHEDRON)
                                continue;

                            for (Integer m = 0; m < anOtherElement.nbNodeIds(); m++) {

                                Integer aNodeId = anOtherElement.nodeId(m);

                                if (aNodeId != aNodeId1 && aNodeId != aNodeId2 && aNodeId != aNodeId3 && aNodeId != aNodeId4) {
                                    aNodeId5 = aNodeId;
                                    bNodeFound = true;
                                    break;
                                }
                            }

                            if (bNodeFound)
                                break;
                        }

                        if (!bNodeFound)
                            continue;

                        CMshNode<Real>& aNode1 = m_volumeMesh.node(aNodeId1);
                        CMshNode<Real>& aNode2 = m_volumeMesh.node(aNodeId2);
                        CMshNode<Real>& aNode3 = m_volumeMesh.node(aNodeId3);
                        CMshNode<Real>& aNode4 = m_volumeMesh.node(aNodeId4);
                        CMshNode<Real>& aNode5 = m_volumeMesh.node(aNodeId5);

                        // Original
                        CMshTetrahedron<Real> aTetrahedron1;
                        aTetrahedron1.addVertex(aNode1);
                        aTetrahedron1.addVertex(aNode2);
                        aTetrahedron1.addVertex(aNode3);
                        aTetrahedron1.addVertex(aNode4);

                        aTetrahedron1.calculateQuality();
                        Real q1 = aTetrahedron1.quality();

                        CMshTetrahedron<Real> aTetrahedron2;
                        aTetrahedron2.addVertex(aNode1);
                        aTetrahedron2.addVertex(aNode5);
                        aTetrahedron2.addVertex(aNode4);
                        aTetrahedron2.addVertex(aNode3);

                        aTetrahedron2.calculateQuality();
                        Real q2 = aTetrahedron2.quality();

                        CMshTetrahedron<Real> aTetrahedron3;
                        aTetrahedron3.addVertex(aNode5);
                        aTetrahedron3.addVertex(aNode2);
                        aTetrahedron3.addVertex(aNode4);
                        aTetrahedron3.addVertex(aNode3);

                        aTetrahedron3.calculateQuality();
                        Real q3 = aTetrahedron3.quality();

                        // Modified
                        CMshTetrahedron<Real> aTetrahedron4;
                        aTetrahedron4.addVertex(aNode1);
                        aTetrahedron4.addVertex(aNode2);
                        aTetrahedron4.addVertex(aNode5);
                        aTetrahedron4.addVertex(aNode4);

                        aTetrahedron4.calculateVolume();
                        aTetrahedron4.calculateQuality();
                        Real q4 = aTetrahedron4.quality();

                        CMshTetrahedron<Real> aTetrahedron5;
                        aTetrahedron5.addVertex(aNode1);
                        aTetrahedron5.addVertex(aNode2);
                        aTetrahedron5.addVertex(aNode3);
                        aTetrahedron5.addVertex(aNode5);

                        aTetrahedron5.calculateVolume();
                        aTetrahedron5.calculateQuality();
                        Real q5 = aTetrahedron5.quality();

                        if (std::min(q4, q5) > std::min(q1, std::min(q2, q3)) && aTetrahedron4.volume() > aTolerance * aTolerance * aTolerance && aTetrahedron5.volume() > aTolerance * aTolerance * aTolerance) {

                            Integer anElementId1 = sElementIds[0];
                            Integer anElementId2 = sElementIds[1];
                            Integer anElementId3 = sElementIds[2];

                            // Do flip
                            m_volumeMesh.element(anElementId1).setNodeId(0, aNodeId1);
                            m_volumeMesh.element(anElementId1).setNodeId(1, aNodeId2);
                            m_volumeMesh.element(anElementId1).setNodeId(2, aNodeId5);
                            m_volumeMesh.element(anElementId1).setNodeId(3, aNodeId4);

                            m_volumeMesh.element(anElementId2).setNodeId(0, aNodeId1);
                            m_volumeMesh.element(anElementId2).setNodeId(1, aNodeId2);
                            m_volumeMesh.element(anElementId2).setNodeId(2, aNodeId3);
                            m_volumeMesh.element(anElementId2).setNodeId(3, aNodeId5);

                            // Contract one of them
                            m_volumeMesh.element(anElementId3).setNodeId(0, aNodeId1);
                            m_volumeMesh.element(anElementId3).setNodeId(1, aNodeId1);
                            m_volumeMesh.element(anElementId3).setNodeId(2, aNodeId1);
                            m_volumeMesh.element(anElementId3).setNodeId(3, aNodeId1);

                            sFlipped[anElementId1] = true;
                            sFlipped[anElementId2] = true;
                            sFlipped[anElementId3] = true;

                            bDeleteElement[anElementId3] = true;
                        }
                    }
                }
            }
        }

        for (typename std::map<Integer, bool>::const_iterator itr = bDeleteElement.begin(); itr != bDeleteElement.end(); ++itr) {

            Integer anElementId = itr->first;

            if (itr->second)
                m_volumeMesh.removeElement(anElementId);
        }

        m_volumeMesh.generateFaces(aTolerance);
    }
}
}
