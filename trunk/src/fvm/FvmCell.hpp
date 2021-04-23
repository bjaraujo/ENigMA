// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <map>

#include <Eigen/Dense>

#include "GeoVector.hpp"

namespace ENigMA
{
    namespace fvm
    {
        template <typename Real>
        class CFvmCell
        {
        public:
            CFvmCell();
            virtual ~CFvmCell();
        };
    }
}

#include "FvmCell_Imp.hpp"
