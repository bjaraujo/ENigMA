// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <map>
#include <vector>

#include "GeoHashGrid.hpp"

#include "FvmMesh.hpp"

namespace ENigMA
{
    namespace fvm
    {
        template <typename Real>
        class CFvmMeshSearch
        {
        private:
            CFvmMesh<Real>* m_mesh;

            CGeoHashGrid<Real> m_boundaryFaceHashGrid;

        public:
            CFvmMeshSearch();
            explicit CFvmMeshSearch(CFvmMesh<Real>& aMesh);
            virtual ~CFvmMeshSearch();

            void set(CFvmMesh<Real>& aMesh);

            void build();

            void findClosestBoundaryFaces(CGeoCoordinate<Real>& aCoordinate, const Real aDistance, std::vector<Integer>& sFaceIds);
        };
    }
}

#include "FvmMeshSearch_Imp.hpp"
