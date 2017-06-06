// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#pragma once

#include "GeoQuadrilateral.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace mesh
    {

        template <typename Real>
        class CMshQuadrilateral : public CGeoQuadrilateral<Real>
        {
        private:

            Real m_quality;

        public:
            CMshQuadrilateral();
            ~CMshQuadrilateral();

            void calculateQuality();
            Real quality();

        };

    }

}

#include "MshQuadrilateral_Imp.hpp"

