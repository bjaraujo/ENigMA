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

        namespace thermal
        {

            template <typename Real>
            class CFemThermalElement
            {
            private:

                Real m_specificHeat;
                Real m_density;
                Real m_thermalConductivity;

            public:

                CFemThermalElement();
                ~CFemThermalElement();

                Real& specificHeat();
                Real& density();
                Real& thermalConductivity();
    
            };

        }

    }

}

#include "FemThermalElement_Imp.hpp"
