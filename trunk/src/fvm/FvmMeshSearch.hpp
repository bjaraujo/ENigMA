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

namespace ENigMA {

namespace fvm {

    template <typename Real>
    class CFvmMeshSearch {
    private:
        CFvmMesh<Real>* m_mesh;

        CGeoHashGrid<Real> m_boundaryFaceHashGrid;

    public:
        CFvmMeshSearch();
        CFvmMeshSearch(CFvmMesh<Real>& aMesh);
        virtual ~CFvmMeshSearch();

        void set(CFvmMesh<Real>& aMesh);

        void build();

        void findClosestBoundaryFace(CGeoCoordinate<Real>& aCoordinate, Integer& aFaceId, const Real aTolerance);
    };
}
}

#include "FvmMeshSearch_Imp.hpp"
