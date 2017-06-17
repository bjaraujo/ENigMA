// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoVolume.hpp"
#include "GeoCoordinate.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoSphere : public CGeoVolume<Real>
        {
        private:

            CGeoCoordinate<Real> m_center;
            Real m_radius;

        public:
            CGeoSphere(const CGeoCoordinate<Real>& aPoint1, const CGeoCoordinate<Real>& aPoint2, const CGeoCoordinate<Real>& aPoint3, const CGeoCoordinate<Real>& aPoint4, const Real aTolerance = 0.0);
            CGeoSphere(const CGeoCoordinate<Real>& aCenter, const Real aRadius);
            ~CGeoSphere();

            CGeoCoordinate<Real>& center();
            Real radius();

            inline void calculateCentroid(bool bReCalculate = false);
            inline void calculateSurfaceArea(bool bReCalculate = false);
            inline void calculateVolume(bool bReCalculate = false);
            inline void calculateBoundingBox(bool bReCalculate = false);

            inline bool contains(const CGeoCoordinate<Real>& aPoint, const Real aTolerance = 0.0);

        };

    }

}

#include "GeoSphere_Imp.hpp"

