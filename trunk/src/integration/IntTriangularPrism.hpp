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
        class CIntTriangularPrism : public CIntGaussIntegration<Real>
        {
        protected:

            std::vector<Real> m_xi, m_eta, m_zeta;
            std::vector<Real> m_wxi, m_weta, m_wzeta;

            void setGaussPoints();

        public:

            CIntTriangularPrism();
            ~CIntTriangularPrism();

        };

    }

}

#include "IntTriangularPrism_Imp.hpp"
