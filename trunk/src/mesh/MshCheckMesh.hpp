// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "MshMesh.hpp"

namespace ENigMA {

namespace mesh {

    template <typename Real>
    class CMshCheckMesh {
    public:
        CMshCheckMesh();
        ~CMshCheckMesh();

        bool checkOpen(CMshMesh<Real>& aMesh);
        bool checkIntersections(CMshMesh<Real>& aMesh, const Real aTolerance = 0.0);
    };
}
}

#include "MshCheckMesh_Imp.hpp"
