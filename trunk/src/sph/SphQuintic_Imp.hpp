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

namespace ENigMA
{

    namespace sph
    {

        template <typename Real>
        CSphQuintic<Real>::CSphQuintic(const Integer nDimension) : CSphKernel<Real>(nDimension)
        {

        }

        template <typename Real>
        CSphQuintic<Real>::CSphQuintic() : CSphKernel<Real>()
        {

        }

        template <typename Real>
        CSphQuintic<Real>::~CSphQuintic()
        {

        }

        template <typename Real>
        void CSphQuintic<Real>::setDimension(const Integer nDimension)
        {

            CSphKernel<Real>::m_dim = nDimension;

        }

        template <typename Real>
        Real CSphQuintic<Real>::W(const CGeoVector<Real> r, const Real h)
        {

            Real q = r.norm() / h;

            if (q > 2.0)
                return 0.0;

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 3.0 / (2.0 * h);
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 7.0 / (4.0 * CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 21.0 / (16.0 * CSphKernel<Real>::m_pi * h * h * h);

            return CSphKernel<Real>::m_C * (pow(1 - 0.5 * q, 4) * (2 * q + 1));

        }

        template <typename Real>
        CGeoVector<Real> CSphQuintic<Real>::gradientW(const CGeoVector<Real> r, const Real h, const Real aTolerance)
        {

            Real q = r.norm() / h;

            if (q < aTolerance)
                q = aTolerance;

            if (q > 2.0)
                return typename CGeoVector<Real>::CGeoVector(0, 0, 0);

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 3.0 / (2.0 * h);
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 7.0 / (4.0 * CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 21.0 / (16.0 * CSphKernel<Real>::m_pi * h * h * h);

            return -CSphKernel<Real>::m_C * (2 * pow(1 - 0.5 * q, 4) - 2 * pow(1 - 0.5 * q, 3) * (2 * q + 1)) / q * r;

        }

        template <typename Real>
        Real CSphQuintic<Real>::laplacianW(const CGeoVector<Real> r, const Real h)
        {

            Real q = r.norm() / h;

            if (q > 2.0)
                return 0.0;

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 3.0 / (2.0 * h);
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 7.0 / (4.0 * CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 21.0 / (16.0 * CSphKernel<Real>::m_pi * h * h * h);

            Real w = pow(1 - 0.5 * q, 4) * (2 * q + 1);
            Real dw = 2 * pow(1 - 0.5 * q, 4) - 2 * pow(1 - 0.5 * q, 3) * (2 * q + 1);

            return CSphKernel<Real>::m_C * h * (-dw * q + w * CSphKernel<Real>::m_dim);

        }

    }

}
