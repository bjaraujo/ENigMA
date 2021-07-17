// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoTriangularPrism.hpp"

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        class CMshTriangularPrism : public ENigMA::geometry::CGeoTetrahedron<Real>
        {
        private:
            Real m_quality;

        public:
            CMshTriangularPrism();
            virtual ~CMshTriangularPrism();

            void calculateQuality();
            Real quality() const;
        };
    }
}

#include "MshTriangularPrism_Imp.hpp"
