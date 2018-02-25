// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FemBeam.hpp"
#include "FemThermalElement.hpp"

namespace ENigMA
{

    namespace fem
    {

        namespace thermal
        {

            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemLinearTemperatureBeam : public CFemThermalElement<Real>, public CFemBeam<Real, NbNodes, Dof, Order>
            {
            };

            template <typename Real>
            class CFemLinearTemperatureBeam<Real, 2, 1, 1> : public CFemThermalElement<Real>, public CFemBeam<Real, 2, 1, 1>
            {
            public:
                
                CFemLinearTemperatureBeam();
                ~CFemLinearTemperatureBeam();

                void setConvectionOnEdge(const Real h, const Real Tinf);
                void setConvectionOnEdge(const Real e, const Real teta, const Real Tinf);

            };

        }

    }

}

#include "FemLinearTemperatureBeam_Imp.hpp"


