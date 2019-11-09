// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "SphKernel.hpp"

namespace ENigMA {

namespace sph {

    template <typename Real>
    class CSphSpiky : public CSphKernel<Real> {
    public:
        CSphSpiky();
        CSphSpiky(const Integer nDimension);
        virtual ~CSphSpiky();

        void setDimension(const Integer nDimension);

        Real W(const ENigMA::geometry::CGeoVector<Real> r, const Real h);
        ENigMA::geometry::CGeoVector<Real> gradientW(const ENigMA::geometry::CGeoVector<Real> r, const Real h, const Real aTolerance = 0.0);
        Real laplacianW(const ENigMA::geometry::CGeoVector<Real> r, const Real h);
    };
}
}

#include "SphSpiky_Imp.hpp"
