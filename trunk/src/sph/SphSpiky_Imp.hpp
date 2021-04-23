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
        CSphSpiky<Real>::CSphSpiky(const Integer nDimension)
            : CSphKernel<Real>(nDimension)
        {
        }

        template <typename Real>
        CSphSpiky<Real>::CSphSpiky()
            : CSphKernel<Real>()
        {
        }

        template <typename Real>
        CSphSpiky<Real>::~CSphSpiky()
        {
        }

        template <typename Real>
        void CSphSpiky<Real>::setDimension(const Integer nDimension)
        {

            CSphKernel<Real>::m_dim = nDimension;
        }

        template <typename Real>
        Real CSphSpiky<Real>::W(const CGeoVector<Real> r, const Real h)
        {

            Real q = r.norm() / h;

            if (q > 1.0)
                return 0.0;

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 4.0 / h;
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 10.0 / (CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 15.0 / (CSphKernel<Real>::m_pi * h * h * h);

            return CSphKernel<Real>::m_C * pow(1 - q, 3);
        }

        template <typename Real>
        CGeoVector<Real> CSphSpiky<Real>::gradientW(const CGeoVector<Real> r, const Real h, const Real aTolerance)
        {

            Real q = r.norm() / h;

            if (q < aTolerance)
                q = aTolerance;

            if (q > 1.0)
                return CGeoVector<Real>(0, 0, 0);

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 4.0 / h;
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 10.0 / (CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 15.0 / (CSphKernel<Real>::m_pi * h * h * h);

            return CSphKernel<Real>::m_C * 3 * pow(1 - q, 2) / q * r;
        }

        template <typename Real>
        Real CSphSpiky<Real>::laplacianW(const CGeoVector<Real> r, const Real h)
        {

            Real q = r.norm() / h;

            if (q > 1.0)
                return 0.0;

            if (CSphKernel<Real>::m_dim == 1)
                CSphKernel<Real>::m_C = 4.0 / h;
            else if (CSphKernel<Real>::m_dim == 2)
                CSphKernel<Real>::m_C = 10.0 / (CSphKernel<Real>::m_pi * h * h);
            else if (CSphKernel<Real>::m_dim == 3)
                CSphKernel<Real>::m_C = 15.0 / (CSphKernel<Real>::m_pi * h * h * h);

            Real w = pow(1 - q, 3);
            Real dw = -3.0 * pow(1 - q, 2);

            return CSphKernel<Real>::m_C * h * (-dw * q + w * CSphKernel<Real>::m_dim);
        }
    }
}
