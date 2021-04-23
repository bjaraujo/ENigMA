// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <cmath>

namespace ENigMA
{
    namespace sph
    {
        template <typename Real>
        CSphKernel<Real>::CSphKernel(const Integer nDimension)
            : m_pi(std::acos(-1.0))
        {
            this->setDimension(nDimension);
        }

        template <typename Real>
        CSphKernel<Real>::CSphKernel()
            : m_pi(std::acos(-1.0))
        {
            this->setDimension(1);
        }

        template <typename Real>
        CSphKernel<Real>::~CSphKernel()
        {
        }

        template <typename Real>
        void CSphKernel<Real>::setDimension(const Integer nDimension)
        {
            m_dim = nDimension;
        }
    }
}
