// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoHexahedron.hpp"
#include "GeoLine.hpp"
#include "GeoQuadrilateral.hpp"
#include "GeoTetrahedron.hpp"
#include "GeoTriangle.hpp"
#include "MshMesh.hpp"

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
            virtual ~CMshBasicMesher();

            bool generate(const ENigMA::geometry::CGeoLine<Real>& aLine, const Integer nu);
            bool generate(const ENigMA::geometry::CGeoQuadrilateral<Real>& aQuadrilateral, const Integer nu, const Integer nv, bool decimate = false);
            bool generate(const ENigMA::geometry::CGeoHexahedron<Real>& aHexahedron, const Integer nu, const Integer nv, const Integer nw, bool decimate = false);

            bool generate(const ENigMA::geometry::CGeoBoundingBox<Real>& aBoundingBox, const Real meshSize, bool decimate = false);

            ENigMA::mesh::CMshMesh<Real>& mesh();
        };
    }
}

#include "MshBasicMesher_Imp.hpp"
