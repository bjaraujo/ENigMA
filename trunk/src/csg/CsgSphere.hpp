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
        template <typename Real>
        class CCsgSphere
        {
        public:
            static ENigMA::geometry::CGeoCoordinate<Real> vertex(ENigMA::geometry::CGeoCoordinate<Real>& aCenter, Real aRadius, Real theta, Real phi);
            static ENigMA::csg::CCsgBoolean<Real> create(ENigMA::geometry::CGeoCoordinate<Real>& aCenter, Real aRadius, Integer nSlices, Integer nStacks);
        };
    }
}

#include "CsgSphere_Imp.hpp"
