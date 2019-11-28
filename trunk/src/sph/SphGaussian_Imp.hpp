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

using namespace ENigMA::geometry;

namespace ENigMA {
namespace sph {
    template <typename Real>
    CSphGaussian<Real>::CSphGaussian(const Integer nDimension)
        : CSphKernel<Real>(nDimension)
    {
    }

    template <typename Real>
    CSphGaussian<Real>::CSphGaussian()
        : CSphKernel<Real>()
    {
    }

    template <typename Real>
    CSphGaussian<Real>::~CSphGaussian()
    {
    }

    template <typename Real>
    void CSphGaussian<Real>::setDimension(const Integer nDimension)
    {
        CSphKernel<Real>::m_dim = nDimension;
    }

    template <typename Real>
    Real CSphGaussian<Real>::W(const CGeoVector<Real> r, const Real h)
    {
        Real q = r.norm() / h;

        if (q > 3.0)
            return 0.0;

        if (CSphKernel<Real>::m_dim == 1)
            CSphKernel<Real>::m_C = 1.0 / (sqrt(CSphKernel<Real>::m_pi) * h);
        else if (CSphKernel<Real>::m_dim == 2)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * h * h);
        else if (CSphKernel<Real>::m_dim == 3)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * CSphKernel<Real>::m_pi * CSphKernel<Real>::m_pi * h * h * h);

        return CSphKernel<Real>::m_C * exp(-q * q);
    }

    template <typename Real>
    CGeoVector<Real> CSphGaussian<Real>::gradientW(const CGeoVector<Real> r, const Real h, const Real aTolerance)
    {
        Real q = r.norm() / h;

        if (q < aTolerance)
            q = aTolerance;

        if (q > 3.0)
            return CGeoVector<Real>(0, 0, 0);

        if (CSphKernel<Real>::m_dim == 1)
            CSphKernel<Real>::m_C = 1.0 / (sqrt(CSphKernel<Real>::m_pi) * h);
        else if (CSphKernel<Real>::m_dim == 2)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * h * h);
        else if (CSphKernel<Real>::m_dim == 3)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * CSphKernel<Real>::m_pi * CSphKernel<Real>::m_pi * h * h * h);

        return CSphKernel<Real>::m_C * 2 * exp(-q * q) * r;
    }

    template <typename Real>
    Real CSphGaussian<Real>::laplacianW(const CGeoVector<Real> r, const Real h)
    {
        Real q = r.norm() / h;

        if (q > 3.0)
            return 0.0;

        if (CSphKernel<Real>::m_dim == 1)
            CSphKernel<Real>::m_C = 1.0 / (sqrt(CSphKernel<Real>::m_pi) * h);
        else if (CSphKernel<Real>::m_dim == 2)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * h * h);
        else if (CSphKernel<Real>::m_dim == 3)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * CSphKernel<Real>::m_pi * CSphKernel<Real>::m_pi * h * h * h);

        Real w = exp(-q * q);
        Real dw = -2.0 * q * w;

        return CSphKernel<Real>::m_C * h * (-dw * q + w * CSphKernel<Real>::m_dim);
    }
}
}
