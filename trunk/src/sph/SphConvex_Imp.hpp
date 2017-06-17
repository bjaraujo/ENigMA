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
        CSphConvex<Real>::CSphConvex(const Integer nDimension) : CSphKernel<Real>(nDimension)
        {

        }

        template <typename Real>
        CSphConvex<Real>::CSphConvex() : CSphKernel<Real>()
        {

        }

        template <typename Real>
        CSphConvex<Real>::~CSphConvex()
        {

        }

        template <typename Real>
        void CSphConvex<Real>::setDimension(const Integer nDimension)
        {

            CSphKernel<Real>::m_dim = nDimension;

        }

        template <typename Real>
        Real CSphConvex<Real>::W(const CGeoVector<Real> r, const Real h)
        {

            Real q = r.norm() / h;

            if (q > 2.0)
                return 0.0;

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 5.0 / (24.0 * h);
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 15.0 / (64.0 * CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 21.0 / (128.0 * CSphKernel<Real>::m_pi * h * h * h);

            return CSphKernel<Real>::m_C * pow(2 - q, 3) * (0.5 * q + 1);

        }

        template <typename Real>
        CGeoVector<Real> CSphConvex<Real>::gradientW(const CGeoVector<Real> r, const Real h, const Real aTolerance)
        {

            Real q = r.norm() / h;

            if (q < aTolerance)
                q = aTolerance;

            if (q > 2.0)
                return typename CGeoVector<Real>::CGeoVector(0, 0, 0);

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 5.0 / (24.0 * h);
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 15.0 / (64.0 * CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 21.0 / (128.0 * CSphKernel<Real>::m_pi * h * h * h);

            return -CSphKernel<Real>::m_C * (pow(2 - q, 1.5) - 3 * (2 - q) * (2 - q) * (0.5 * q + 1)) / q * r;

        }

        template <typename Real>
        Real CSphConvex<Real>::laplacianW(const CGeoVector<Real> r, const Real h)
        {

            Real q = r.norm() / h;

            if (q > 2.0)
                return 0.0;

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 5.0 / (24.0 * h);
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 15.0 / (64.0 * CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 21.0 / (128.0 * CSphKernel<Real>::m_pi * h * h * h);

            Real w = pow(2 - q, 3) * (0.5 * q + 1);
            Real dw = pow(2 - q, 1.5) - 3 * (2 - q) * (2 - q) * (0.5 * q + 1);

            return CSphKernel<Real>::m_C * h * (dw * q + w * CSphKernel<Real>::m_dim);

        }

    }

}
