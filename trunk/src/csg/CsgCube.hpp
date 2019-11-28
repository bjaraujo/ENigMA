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
        class CCsgCube
        {
        public:

            static ENigMA::csg::CCsgBoolean<Real> create(ENigMA::geometry::CGeoCoordinate<Real> aCenter, ENigMA::geometry::CGeoVector<Real> aRadius);

        };
    }         
}
    
#include "CsgCube_Imp.hpp"
