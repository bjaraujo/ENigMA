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
    struct SMshQuadrilateralAdvancingFrontEdge {
        Integer id;

        bool remove;
        bool boundary;

        Integer nodeId[2];

        Integer neighborId[2];

        Integer quadrilateralId;

        CGeoLine<Real> line;

        void build(CMshMesh<Real>& aMesh)
        {
            this->line.reset();
            this->line.setStartPoint(aMesh.node(this->nodeId[0]));
            this->line.setEndPoint(aMesh.node(this->nodeId[1]));
        }
    };

    template <typename Real>
    class CMshQuadrilateralMesher {
    private:
        // Inner points
        struct SNode {
            Integer id;
            bool remove;

            Integer nodeId;
        };

        std::vector<SNode> m_interiorNodes;

        std::vector<SMshQuadrilateralAdvancingFrontEdge<Real>> m_anAdvFront;

        CGeoBoundingBox<Real> m_boundingBox;

        clock_t m_begin;
        clock_t m_end;

        bool m_bStop;
        Integer m_dataInterval;
        Integer m_timeInterval;
        Integer m_previousNbElements;

        Integer m_nextEdgeId;

        CGeoRtree<Real> m_tree;
        CMshMesh<Real> m_surfaceMesh;

        void checkUpdate();

        inline bool edgeExists(SMshQuadrilateralAdvancingFrontEdge<Real>& anAdvEdge, Integer& aDuplicateEdgeId, std::vector<Integer>& sEdges, const Real aTolerance = 0.0);
        inline bool edgeOk(SMshQuadrilateralAdvancingFrontEdge<Real>& anAdvEdge, CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, std::vector<Integer>& sEdges, const Real aTolerance = 0.0);
        inline bool quadrilateralContainsNode(CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, CMshNode<Real>& aNode4, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance = 0.0);
        inline bool checkDelaunay(CMshNode<Real>& aNewNode, const Real aTolerance = 0.0);

        void cleanDuplicateEdges(std::vector<Integer>& sEdges, const Real aTolerance = 0.0);
        void adjustConnectivity(std::vector<Integer>& sEdges);
        void findClosestNodes(std::vector<Integer>& sEdges, std::vector<Integer>& sNodes);
        Real findShortestDistance(std::vector<Integer>& sEdges, CGeoLine<Real>& aLine, Integer anAdvEdgeId, const Real aTolerance);

        void addQuadrilateral(SMshQuadrilateralAdvancingFrontEdge<Real>& anAdvEdge, const Integer aNodeId3, const Integer aNodeId4, std::vector<Integer>& sEdges, const Real aTolerance);

        void addEdgeToRtree(SMshQuadrilateralAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);
        void removeEdgeFromRtree(SMshQuadrilateralAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);

        void removeEdge(SMshQuadrilateralAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);

        bool advancingFrontQuadMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real sizeFactor = 1.0, Real shrinkFactor = 1.0, Real expandFactor = 1.0, Real minQuality = 0.0, const bool bCheckDelaunay = false, Integer firstIndex = 0, const Real aTolerance = 0.0);

        Integer frontSize();

    public:
        CMshQuadrilateralMesher();
        virtual ~CMshQuadrilateralMesher();

        bool remesh(CMshMesh<Real>& anEdgeMesh, Real meshSize);
        bool remesh(CMshMesh<Real>& anEdgeMesh, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc);

        bool generate(CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, Real meshSize, Real minQuality = 0.0, const Real aTolerance = 0.0);
        bool generate(CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, Real meshSize, Real minQuality = 0.0, const Real aTolerance = 0.0);
        bool generate(CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minQuality = 0.0, const Real aTolerance = 0.0);

        CMshMesh<Real>& mesh();

        void setIntervals(const Integer timeInterval, const Integer dataInterval);
        void stopMeshing();

        void applyFixedBoundary(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Real aTolerance = 0.0);

        void flipEdges(const Real aTolerance = 0.0);
        void relaxNodes(const Real aTolerance = 0.0);

        // callback
        std::function<int(int)> onUpdate;
    };
}
}

#include "MshQuadrilateralMesher_Imp.hpp"
