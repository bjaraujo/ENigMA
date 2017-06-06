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
        class CIntTriangle : public CIntGaussIntegration<Real>
        {
        protected:

            std::vector<Real> m_xi, m_eta;
            std::vector<Real> m_wxi, m_weta;

            std::vector<Real> m_beta1, m_beta2, m_beta3;    // barycentric coordinates
            std::vector<Real> m_wbeta;                        // weights

            void setGaussPoints();
            void setBarycentricGaussPoints();

        public:

            CIntTriangle();
            ~CIntTriangle();

        };

    }

}

#include "IntTriangle_Imp.hpp"
