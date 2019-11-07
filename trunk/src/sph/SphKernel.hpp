// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "CmnTypes.hpp"
#include "GeoVector.hpp"

namespace ENigMA {

namespace sph {

    template <typename Real>
    class CSphKernel {
    protected:
        Integer m_dim;
        Real m_pi;
        Real m_C;

    public:
        CSphKernel();
        CSphKernel(const Integer nDimension);
        ~CSphKernel();

        virtual void setDimension(const Integer nDimension);

        virtual Real W(const ENigMA::geometry::CGeoVector<Real> r, const Real h) = 0;
        virtual ENigMA::geometry::CGeoVector<Real> gradientW(const ENigMA::geometry::CGeoVector<Real> r, const Real h, const Real aTolerance = 0.0) = 0;
        virtual Real laplacianW(const ENigMA::geometry::CGeoVector<Real> r, const Real h) = 0;
    };
}
}

#include "SphKernel_Imp.hpp"
