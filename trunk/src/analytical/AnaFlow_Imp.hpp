// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {

namespace analytical {

    template <typename Real>
    CAnaFlow<Real>::CAnaFlow()
    {
    }

    template <typename Real>
    CAnaFlow<Real>::~CAnaFlow()
    {
    }

    template <typename Real>
    Real CAnaFlow<Real>::phi(const Real x, const Real t, const Real L, const Real mu)
    {

        return exp(-(x - 4 * t) * (x - 4 * t) / (4 * mu * (t + 1))) + exp(-(x - 4 * t - L) * (x - 4 * t - L) / (4 * mu * (t + 1)));
    }

    template <typename Real>
    void CAnaFlow<Real>::viscousBurgersEquation(Real x, Real dx, Real t, Real L, Real mu, Real& u)
    {

        Real phi0 = phi(x, t, L, mu);
        Real phi1 = phi(x + dx, t, L, mu);

        u = -2 * mu / phi0 * (phi1 - phi0) / dx + 4;
    }
}
}
