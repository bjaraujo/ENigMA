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

        namespace structural
        {

            template <typename Real>
            class CFemStructuralElement
            {
            private:

                Real m_coeffPoisson;
                Real m_elasticModulus;

            public:

                CFemStructuralElement();
                ~CFemStructuralElement();

                Real& coeffPoisson();
                Real& elasticModulus();

            };

        }

    }

}

#include "FemStructuralElement_Imp.hpp"
