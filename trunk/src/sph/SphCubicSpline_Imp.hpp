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
    CSphCubicSpline<Real>::CSphCubicSpline(const Integer nDimension)
        : CSphKernel<Real>(nDimension)
    {
    }

    template <typename Real>
    CSphCubicSpline<Real>::CSphCubicSpline()
        : CSphKernel<Real>()
    {
    }

    template <typename Real>
    CSphCubicSpline<Real>::~CSphCubicSpline()
    {
    }

    template <typename Real>
    void CSphCubicSpline<Real>::setDimension(const Integer nDimension)
    {

        CSphKernel<Real>::m_dim = nDimension;
    }

    template <typename Real>
    Real CSphCubicSpline<Real>::W(const CGeoVector<Real> r, const Real h)
    {

        Real q = r.norm() / h;

        if (q > 2.0)
            return 0.0;

        if (CSphKernel<Real>::m_dim == 1)
            CSphKernel<Real>::m_C = 2.0 / (3.0 * h);
        else if (CSphKernel<Real>::m_dim == 2)
            CSphKernel<Real>::m_C = 10.0 / (7.0 * CSphKernel<Real>::m_pi * h * h);
        else if (CSphKernel<Real>::m_dim == 3)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * h * h * h);

        if (q >= 0 && q < 1.0)
            return CSphKernel<Real>::m_C * ((2 - q) * (2 - q) * (2 - q) - 4 * (1 - q) * (1 - q) * (1 - q));
        else
            return CSphKernel<Real>::m_C * (2 - q) * (2 - q) * (2 - q);
    }

    template <typename Real>
    CGeoVector<Real> CSphCubicSpline<Real>::gradientW(const CGeoVector<Real> r, const Real h, const Real aTolerance)
    {

        Real q = r.norm() / h;

        if (q < aTolerance)
            q = aTolerance;

        if (q > 2.0)
            return typename CGeoVector<Real>::CGeoVector(0, 0, 0);

        if (CSphKernel<Real>::m_dim == 1)
            CSphKernel<Real>::m_C = 2.0 / (3.0 * h);
        else if (CSphKernel<Real>::m_dim == 2)
            CSphKernel<Real>::m_C = 10.0 / (7.0 * CSphKernel<Real>::m_pi * h * h);
        else if (CSphKernel<Real>::m_dim == 3)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * h * h * h);

        if (q >= aTolerance && q < 1.0)
            return CSphKernel<Real>::m_C * (3.0 * q * (1 - 0.75 * q)) / q * r;
        else
            return CSphKernel<Real>::m_C * 0.75 * (2 - q) * (2 - q) / q * r;
    }

    template <typename Real>
    Real CSphCubicSpline<Real>::laplacianW(const CGeoVector<Real> r, const Real h)
    {

        Real q = r.norm() / h;

        if (q > 2.0)
            return 0.0;

        if (CSphKernel<Real>::m_dim == 1)
            CSphKernel<Real>::m_C = 2.0 / (3.0 * h);
        else if (CSphKernel<Real>::m_dim == 2)
            CSphKernel<Real>::m_C = 10.0 / (7.0 * CSphKernel<Real>::m_pi * h * h);
        else if (CSphKernel<Real>::m_dim == 3)
            CSphKernel<Real>::m_C = 1.0 / (CSphKernel<Real>::m_pi * h * h * h);

        Real w, dw;

        if (q >= 0 && q < 1.0) {
            w = 1.0 - 1.5 * q * q * (1.0 - 0.5 * q);
            dw = -3.0 * q * (1.0 - 0.75 * q);
        } else {
            w = 0.25 * (2 - q) * (2 - q) * (2 - q);
            dw = -0.75 * (2 - q) * (2 - q);
        }

        return CSphKernel<Real>::m_C * h * (-dw * q + w * CSphKernel<Real>::m_dim);
    }
}
}
