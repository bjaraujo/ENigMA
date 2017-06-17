// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoTriangle.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace mesh
    {

        template <typename Real>
        class CMshTriangle : public CGeoTriangle<Real>
        {
        private:

            Real m_quality;

        public:
            CMshTriangle();
            ~CMshTriangle();

            void calculateQuality();
            Real quality();

        };

    }

}

#include "MshTriangle_Imp.hpp"

