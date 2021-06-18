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

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        class CGeoPlane
        {
        private:
            CGeoNormal<Real> m_normal;
            Real m_d;

        public:
            CGeoPlane();
            CGeoPlane(const CGeoNormal<Real>& aNormal, const Real d);
            virtual ~CGeoPlane();

            CGeoNormal<Real> normal() const;

            void setD(const Real aValue);
            Real d() const;

            void set(const CGeoNormal<Real>& aNormal, const Real d);

            void invert();

            inline Real distance(const CGeoCoordinate<Real>& aPoint);
            inline Real distance(const CGeoCoordinate<Real>& aPoint, CGeoCoordinate<Real>& aNewPoint);
        };
    }
}

#include "GeoPlane_Imp.hpp"
