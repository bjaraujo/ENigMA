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

namespace ENigMA
{

    namespace fem
    {

        namespace flow
        {

            template <typename Real>
            class CFemFlowElement
            {
            private:

                Real m_density;
                Real m_viscosity;

            public:

                CFemFlowElement();
                ~CFemFlowElement();

                Real& density();
                Real& viscosity();

            };

        }

    }

}

#include "FemFlowElement_Imp.hpp"
