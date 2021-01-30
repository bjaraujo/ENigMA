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
#include "GeoLine.hpp"
#include "GeoRtree.hpp"
#include "MshMesh.hpp"

namespace ENigMA {
namespace mesh {

    // Advancing front
    template <typename Real>
    struct SMshTriangleAdvancingFrontEdge {
        Integer id;

        bool remove;
        bool boundary;

        Integer nodeId[2];

        Integer neighborId[2];

        Integer triangleId;

        Integer nodeNotId3;

        ENigMA::geometry::CGeoLine<Real> line;

        void build(CMshMesh<Real>& aMesh)
        {
            this->line.reset();
            this->line.setStartPoint(aMesh.node(this->nodeId[0]));
            this->line.setEndPoint(aMesh.node(this->nodeId[1]));
        }
    };

    template <typename Real>
    class CMshTriangleMesher {
    private:
        // Inner points
        struct SNode {
            Integer id;
            bool remove;

            Integer nodeId;
        };

        std::vector<SNode> m_interiorNodes;

        std::vector<SMshTriangleAdvancingFrontEdge<Real>> m_anAdvFront;

        ENigMA::geometry::CGeoBoundingBox<Real> m_boundingBox;

        clock_t m_begin;
        clock_t m_end;

        bool m_bStop;
        Integer m_dataInterval;
        Integer m_timeInterval;
        Integer m_previousNbElements;

        Integer m_nextEdgeId;

        ENigMA::geometry::CGeoRtree<Real> m_tree;
        ENigMA::mesh::CMshMesh<Real> m_surfaceMesh;

        void checkUpdate();

        inline bool edgeExists(SMshTriangleAdvancingFrontEdge<Real>& anAdvEdge, Integer& aDuplicateEdgeId, std::vector<Integer>& sEdges);
        inline bool edgeOk(SMshTriangleAdvancingFrontEdge<Real>& anAdvEdge, ENigMA::mesh::CMshNode<Real>& aNode1, ENigMA::mesh::CMshNode<Real>& aNode2, std::vector<Integer>& sEdges, const Real aTolerance = 0.0);
        inline bool triangleContainsNode(ENigMA::mesh::CMshNode<Real>& aNode1, ENigMA::mesh::CMshNode<Real>& aNode2, ENigMA::mesh::CMshNode<Real>& aNode3, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance = 0.0);
        inline bool checkDelaunay(ENigMA::mesh::CMshNode<Real>& aNewNode, const Real aTolerance = 0.0);

        void cleanDuplicateEdges(std::vector<Integer>& sEdges, const Real aTolerance = 0.0);
        void adjustConnectivity(std::vector<Integer>& sEdges);
        void findClosestNodes(std::vector<Integer>& sEdges, std::vector<Integer>& sNodes);
        Real findShortestDistance(std::vector<Integer>& sEdges, ENigMA::geometry::CGeoLine<Real>& aLine, Integer anAdvEdgeId, const Real aTolerance);

        void addTriangle(SMshTriangleAdvancingFrontEdge<Real>& anAdvEdge, const Integer aNodeId, std::vector<Integer>& sEdges, const Real aTolerance = 0.0);

        void addEdgeToRtree(SMshTriangleAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);
        void removeEdgeFromRtree(SMshTriangleAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);

        void removeEdge(SMshTriangleAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);

        bool advancingFrontTriMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), Real sizeFactor = 1.0, Real shrinkFactor = 1.0, Real expandFactor = 1.0, Real minQuality = 0.0, const bool bCheckDelaunay = false, Integer firstIndex = 0, const Real aTolerance = 0.0);

        Integer frontSize();

    public:
        CMshTriangleMesher();
        virtual ~CMshTriangleMesher();

        bool remesh(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, Real meshSize);
        bool remesh(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc);

        bool generate(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, Real meshSize, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0);
        bool generate(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<ENigMA::geometry::CGeoCoordinate<Real>>& sInteriorPoints, Real meshSize, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0);
        bool generate(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<ENigMA::geometry::CGeoCoordinate<Real>>& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0);

        ENigMA::mesh::CMshMesh<Real>& mesh();

        void setIntervals(const Integer timeInterval, const Integer dataInterval);
        void stopMeshing();

        void applyFixedBoundary(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Real aTolerance = 0.0);

        void flipEdges(const Real aTolerance = 0.0);
        void relaxNodes(const Real aTolerance = 0.0);
        void collapseEdges(Real collapseSize, const Real aTolerance = 0.0);

        // callback
        std::function<int(int)> onUpdate;
    };
}
}

#include "MshTriangleMesher_Imp.hpp"
