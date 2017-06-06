// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
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
            ~CFvmCell();
            
        };

    }

}

#include "FvmCell_Imp.hpp"
