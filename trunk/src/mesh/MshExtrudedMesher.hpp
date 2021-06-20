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

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        class CMshExtrudedMesher
        {
        private:
            ENigMA::mesh::CMshMesh<Real> m_mesh;

        public:
            CMshExtrudedMesher();
            virtual ~CMshExtrudedMesher();

            bool generate(CMshMesh<Real>& aSurfaceMesh, const Real aDelta, const Real aTolerance = 0.0);

            ENigMA::mesh::CMshMesh<Real>& mesh();
        };
    }
}

#include "MshExtrudedMesher_Imp.hpp"
