// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoArea.hpp"
#include "GeoCoordinate.hpp"

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        class CGeoCircle : public CGeoArea<Real>
        {
        private:
            CGeoCoordinate<Real> m_center;
            Real m_radius;

        public:
            CGeoCircle(const CGeoCoordinate<Real>& aPoint1, const CGeoCoordinate<Real>& aPoint2, const CGeoCoordinate<Real>& aPoint3, const Real aTolerance = 0.0);
            CGeoCircle(const CGeoCoordinate<Real>& aCenter, const Real aRadius);
            virtual ~CGeoCircle();

            const CGeoCoordinate<Real>& center() const;
            Real radius() const;

            inline void calculateBoundingBox(bool bReCalculate = false);
            inline void calculateCentroid(bool bReCalculate = false);
            inline void calculateNormal(bool bReCalculate = false);
            inline void calculateArea(bool bReCalculate = false);

            inline bool contains(const CGeoCoordinate<Real>& aPoint, const Real aTolerance = 0.0);
        };
    }
}

#include "GeoCircle_Imp.hpp"
