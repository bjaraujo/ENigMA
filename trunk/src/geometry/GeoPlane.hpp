// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoCoordinate.hpp"
#include "GeoNormal.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    class CGeoPlane {
    private:
        CGeoNormal<Real> m_normal;
        Real m_d;

    public:
        CGeoPlane();
        CGeoPlane(CGeoNormal<Real>& aNormal, const Real d);
        virtual ~CGeoPlane();

        CGeoNormal<Real>& normal();

        void setD(const Real aValue);
        Real d() const;

        void set(CGeoNormal<Real>& aNormal, const Real d);

        void invert();

        inline bool distance(const CGeoCoordinate<Real>& aPoint, CGeoCoordinate<Real>& aNewPoint, Real& aDistance, const Real aTolerance = 0.0);
    };
}
}

#include "GeoPlane_Imp.hpp"
