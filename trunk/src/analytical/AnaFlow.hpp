// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "AnaFunction.hpp"

namespace ENigMA {
namespace analytical {
    template <typename Real>
    class CAnaFlow : public CAnaFunction<Real> {
    private:
        Real phi(const Real x, const Real t, const Real L, const Real mu);

    public:
        CAnaFlow();
        virtual ~CAnaFlow();

        void viscousBurgersEquation(Real x, Real dx, Real t, Real L, Real mu, Real& u);
    };
}
}

#include "AnaFlow_Imp.hpp"
