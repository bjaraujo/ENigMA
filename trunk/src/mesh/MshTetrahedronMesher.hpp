// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <ctime>
#include <functional>

#include "AnaFunction.hpp"
#include "GeoRtree.hpp"
#include "GeoTriangle.hpp"
#include "MshMesh.hpp"

namespace ENigMA
{
    namespace mesh
    {

        // Advancing front
        template <typename Real>
        struct SMshAdvancingFrontTriangle
        {
            Integer id;

            bool remove;
            bool boundary;

            Integer nodeId[3];

            Integer neighborId[3];
            Integer nodeNotId[3];

            Integer elementId;
            Integer nodeNotId4;

            CGeoTriangle<Real> triangle;

            void build(const CMshMesh<Real>& aMesh)
            {
                this->triangle.reset();
                this->triangle.addVertex(aMesh.node(this->nodeId[0]));
                this->triangle.addVertex(aMesh.node(this->nodeId[1]));
                this->triangle.addVertex(aMesh.node(this->nodeId[2]));
                this->triangle.calculateArea(true);
            }
        };

        template <typename Real>
        struct SNode
        {
            Integer id;
            bool remove;

            Integer nodeId;
        };

        template <typename Real>
        class CMshTetrahedronMesher
        {
        private:
            std::vector<SNode<Real>> m_innerNodes;

            std::vector<SMshAdvancingFrontTriangle<Real>> m_anAdvFront;

            CGeoBoundingBox<Real> m_boundingBox;

            clock_t m_begin;
            clock_t m_end;

            bool m_bStop;
            Integer m_dataInterval;
            Integer m_timeInterval;
            Integer m_previousNbElements;

            Integer m_nextTriangleId;

            CGeoRtree<Real> m_tree;
            CMshMesh<Real> m_volumeMesh;

            void checkUpdate();

            inline bool triangleExists(SMshAdvancingFrontTriangle<Real>& anAdvTriangle, Integer& aDuplicateTriangleId, std::vector<Integer>& sTriangles);
            inline bool triangleOk(SMshAdvancingFrontTriangle<Real>& anAdvTriangle, CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, std::vector<Integer>& sTriangles, const Real aTolerance = 0.0);
            inline bool tetrahedronContainsNode(CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, CMshNode<Real>& aNode4, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance = 0.0);
            inline bool checkDelaunay(CMshNode<Real>& aNewNode, const Real aTolerance = 0.0);

            bool pairEdges(SMshAdvancingFrontTriangle<Real>& anAdvTriangle1, SMshAdvancingFrontTriangle<Real>& anAdvTriangle2);

            void cleanDuplicateTriangles(std::vector<Integer>& sTriangles, const Real aTolerance = 0.0);
            void adjustConnectivity(std::vector<Integer>& sTriangles);
            void findClosestNodes(std::vector<Integer>& sTriangles, std::vector<Integer>& sNodes);
            Real findShortestDistance(std::vector<Integer>& sTriangles, CGeoLine<Real>& aLine, Integer anAdvTriangleId, const Real aTolerance = 0.0);

            void addTriangleToRtree(SMshAdvancingFrontTriangle<Real>& anAdvTriangle, const Real aTolerance = 0.0);
            void removeTriangleFromRtree(SMshAdvancingFrontTriangle<Real>& anAdvTriangle, const Real aTolerance = 0.0);

            void removeTriangle(SMshAdvancingFrontTriangle<Real>& anAdvTriangle, const Real aTolerance);

            void addTetrahedron(SMshAdvancingFrontTriangle<Real>& anAdvTriangle, const Integer aNodeId, std::vector<Integer>& sTriangles, const Real aTolerance = 0.0);

            bool advancingFrontTetraMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), Real sizeFactor = 1.0, Real shrinkFactor = 1.0, Real expandFactor = 1.0, Real minQuality = 0.0, const bool bAddNodes = true, const bool bCheckDelaunay = false, const Real aTolerance = 0.0);

            Integer getFirstIndex();

            bool repair(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, const Real sizeFactor = 1.0, const Real minQuality = 0.0, const Real aTolerance = 0.0);

            bool rebuildConnectivity(const Real aTolerance = 0.0);

            void reduceFront(const Real aTolerance = 0.0);

            Integer frontSize();

        public:
            CMshTetrahedronMesher();
            virtual ~CMshTetrahedronMesher();

            bool generate(const CMshMesh<Real>& aSurfaceMesh, const Integer maxNbElements, Real meshSize, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0);
            bool generate(const CMshMesh<Real>& aSurfaceMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, Real meshSize, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0);
            bool generate(const CMshMesh<Real>& aSurfaceMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0);

            CMshMesh<Real>& mesh();

            void setIntervals(const Integer timeInterval, const Integer dataInterval);
            void stopMeshing();

            void relaxNodes(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aFactor, const Real aTolerance = 0.0);
            void flipEdges23(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aTolerance = 0.0);
            void flipEdges32(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aTolerance = 0.0);

            // callback
            std::function<bool(bool)> onUpdate;
        };
    }
}

#include "MshTetrahedronMesher_Imp.hpp"
