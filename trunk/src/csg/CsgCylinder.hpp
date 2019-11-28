// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

// Based on Evan Wallace CSG code.

#pragma once

#include "CsgBoolean.hpp"

namespace ENigMA
{
    namespace csg
    {
        template<typename Real> 
        class CCsgCylinder
        {
        public:

            static ENigMA::geometry::CGeoCoordinate<Real> point(ENigMA::geometry::CGeoCoordinate<Real>& aStart, ENigMA::geometry::CGeoVector<Real>& axisX, ENigMA::geometry::CGeoVector<Real>& axisY, ENigMA::geometry::CGeoVector<Real>& axisZ, ENigMA::geometry::CGeoVector<Real>& aRay, Real aRadius, Real stack, Real slice, Real normalBlend);
            static ENigMA::csg::CCsgBoolean<Real> create(ENigMA::geometry::CGeoCoordinate<Real>& aStart, ENigMA::geometry::CGeoCoordinate<Real>& anEnd, Real aRadius, Integer nSlices);

        };
    }
}
    
#include "CsgCylinder_Imp.hpp"
