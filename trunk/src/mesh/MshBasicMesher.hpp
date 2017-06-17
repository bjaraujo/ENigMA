// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoLine.hpp"
#include "GeoTriangle.hpp"
#include "GeoQuadrilateral.hpp"
#include "GeoTetrahedron.hpp"
#include "GeoHexahedron.hpp"
#include "MshMesh.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace mesh
    {

        template <typename Real>
        class CMshBasicMesher
        {
        private:

            ENigMA::mesh::CMshMesh<Real> m_mesh;

        public:

            CMshBasicMesher();
            ~CMshBasicMesher();

            bool generate(CGeoLine<Real>& aLine, const Integer nu);
            bool generate(CGeoQuadrilateral<Real>& aQuadrilateral, const Integer nu, const Integer nv, bool decimate = false);
            bool generate(CGeoHexahedron<Real>& aHexahedron, const Integer nu, const Integer nv, const Integer nw, bool decimate = false);

            ENigMA::mesh::CMshMesh<Real>& mesh();

        };

    }

}

#include "MshBasicMesher_Imp.hpp"
