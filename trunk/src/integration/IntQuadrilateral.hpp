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

#include "IntGaussIntegration.hpp"

namespace ENigMA
{

    namespace integration
    {

        template <typename Real>
        class CIntQuadrilateral : public CIntGaussIntegration<Real>
        {
        protected:

            std::vector<Real> m_xi, m_eta;
            std::vector<Real> m_wxi, m_weta;

            void setGaussPoints();

        public:

            CIntQuadrilateral();
            ~CIntQuadrilateral();

        };

    }

}

#include "IntQuadrilateral_Imp.hpp"
