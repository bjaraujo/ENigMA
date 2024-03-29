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
#include "MshSNode.hpp"

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        struct SMshAdvancingFrontEdge
        {
            Integer id;

            bool remove;
            bool boundary;

            Integer nodeId[2];

            Integer neighborId[2];

            Integer elementId;

            Integer nodeNotId3;

            ENigMA::geometry::CGeoLine<Real> line;

            void build(const CMshMesh<Real>& aMesh)
            {
                this->line.reset();
                this->line.setStartPoint(aMesh.node(this->nodeId[0]));
                this->line.setEndPoint(aMesh.node(this->nodeId[1]));
                this->line.calculateLength();
            }
        };

        template <typename Real>
        class CMshTriangleMesher
        {
        protected:
            std::vector<SNode<Real>> m_interiorNodes;

            std::vector<SMshAdvancingFrontEdge<Real>> m_anAdvFront;
            std::vector<ENigMA::geometry::CGeoLine<Real>> m_anAdvFrontLines;

            ENigMA::geometry::CGeoBoundingBox<Real> m_boundingBox;

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

            inline bool edgeExists(SMshAdvancingFrontEdge<Real>& anAdvEdge, Integer& aDuplicateEdgeId, std::vector<Integer>& sEdges);
            inline bool edgeOk(SMshAdvancingFrontEdge<Real>& anAdvEdge, ENigMA::mesh::CMshNode<Real>& aNode1, ENigMA::mesh::CMshNode<Real>& aNode2, std::vector<Integer>& sEdges, const Real aTolerance = 0.0);
            inline bool triangleContainsNode(ENigMA::mesh::CMshNode<Real>& aNode1, ENigMA::mesh::CMshNode<Real>& aNode2, ENigMA::mesh::CMshNode<Real>& aNode3, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance = 0.0);
            inline bool checkDelaunay(ENigMA::mesh::CMshNode<Real>& aNewNode, const Real aTolerance = 0.0);

            virtual void cleanDuplicateEdges(std::vector<Integer>& sEdges, const Real aTolerance = 0.0);
            void adjustConnectivity(std::vector<Integer>& sEdges);
            void findClosestNodes(std::vector<Integer>& sEdges, std::vector<Integer>& sNodes);
            Real findShortestDistance(std::vector<Integer>& sEdges, ENigMA::geometry::CGeoLine<Real>& aLine, Integer anAdvEdgeId, const Real aTolerance);

            void addTriangle(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Integer aNodeId, std::vector<Integer>& sEdges, const Real aTolerance = 0.0);

            void addEdgeToRtree(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);
            void removeEdgeFromRtree(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);

            void removeEdge(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Real aTolerance = 0.0);

            bool advancingFrontTriMeshing(const Real meshSize, Integer& maxNbElements, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), Real sizeFactor = 1.0, Real shrinkFactor = 1.0, Real expandFactor = 1.0, Real minQuality = 0.0, const bool bAddNodes = true, const bool bCheckDelaunay = false, const Real aTolerance = 0.0);

            Integer frontSize();

        public:
            CMshTriangleMesher();
            virtual ~CMshTriangleMesher();

            virtual bool remesh(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real meshSize, const Real aTolerance = 0.0);

            virtual bool generate(const ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<ENigMA::geometry::CGeoCoordinate<Real>>& sInteriorPoints, const Real meshSize, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0);

            ENigMA::mesh::CMshMesh<Real>& mesh();

            void setIntervals(const Integer timeInterval, const Integer dataInterval);
            void stopMeshing();

            virtual void applyFixedBoundary(ENigMA::mesh::CMshMesh<Real>& aSurfaceMesh, ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Real aTolerance = 0.0);

            virtual void flipEdges(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aTolerance = 0.0);
            virtual void relaxNodes(ENigMA::mesh::CMshMesh<Real>& aMesh, const Real aTolerance = 0.0);
            virtual void collapseEdges(ENigMA::mesh::CMshMesh<Real>& aMesh, Real collapseSize, const Real aTolerance = 0.0);

            // callback
            std::function<int(int)> onUpdate;
        };
    }
}

#include "MshTriangleMesher_Imp.hpp"
