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
#include "MshTriangleMesher.hpp"

namespace ENigMA
{
    namespace mesh
    {

        template <typename Real>
        class CMshQuadrilateralMesher : public CMshTriangleMesher<Real>
        {
        private:
            void cleanDuplicateEdges(std::vector<Integer>& sEdges, const Real aTolerance = 0.0) override;

            inline bool quadrilateralContainsNode(CMshNode<Real>& aNode1, CMshNode<Real>& aNode2, CMshNode<Real>& aNode3, CMshNode<Real>& aNode4, Integer& aNodeId, std::vector<Integer>& sNodes, const Real aTolerance = 0.0);
            void addQuadrilateral(SMshAdvancingFrontEdge<Real>& anAdvEdge, const Integer aNodeId3, const Integer aNodeId4, std::vector<Integer>& sEdges, const Real aTolerance);

            bool advancingFrontQuadMeshing(ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Integer& maxNbElements, Real minMeshSize, Real maxMeshSize, Real sizeFactor = 1.0, Real shrinkFactor = 1.0, Real expandFactor = 1.0, Real minQuality = 0.0, const bool bAddNodes = true, const bool bCheckDelaunay = false, const Real aTolerance = 0.0);

        public:
            CMshQuadrilateralMesher();
            virtual ~CMshQuadrilateralMesher();

            bool remesh(CMshMesh<Real>& anEdgeMesh, Real meshSize) override;
            bool remesh(CMshMesh<Real>& anEdgeMesh, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc) override;

            bool generate(const CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, Real meshSize, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0) override;
            bool generate(const CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, Real meshSize, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0) override;
            bool generate(const CMshMesh<Real>& anEdgeMesh, const Integer maxNbElements, std::vector<CGeoCoordinate<Real>>& sInteriorPoints, ENigMA::analytical::CAnaFunction<Real>& meshSizeFunc, Real minMeshSize = 0.0, Real maxMeshSize = std::numeric_limits<Real>::max(), const Real aTolerance = 0.0) override;

            void applyFixedBoundary(ENigMA::mesh::CMshMesh<Real>& anEdgeMesh, const Real aTolerance = 0.0) override;

            void flipEdges(const Real aTolerance = 0.0) override;
            void relaxNodes(const Real aTolerance = 0.0) override;
        };
    }
}

#include "MshQuadrilateralMesher_Imp.hpp"
