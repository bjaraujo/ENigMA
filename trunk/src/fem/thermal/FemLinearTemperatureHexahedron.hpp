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

#include "../../src/fem/FemHexahedron.hpp"
#include "FemThermalElement.hpp"

namespace ENigMA
{

    namespace fem
    {

        namespace thermal
        {

            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemLinearTemperatureHexahedron : public CFemThermalElement<Real>, public CFemHexahedron<Real, NbNodes, Dof, Order>
            {
            };

            template <typename Real>
            class CFemLinearTemperatureHexahedron<Real, 8, 1, 1> : public CFemThermalElement<Real>, public CFemHexahedron<Real, 8, 1, 1>
            {
            public:
                
                CFemLinearTemperatureHexahedron();
                ~CFemLinearTemperatureHexahedron();

                void setConvectionOnFace(const Integer aFaceIndex, const Real h, const Real Tinf);
                void setConvectionOnFace(const Integer aFaceIndex, const Real e, const Real teta, const Real Tinf);

            };

        }

    }

}

#include "FemLinearTemperatureHexahedron_Imp.hpp"


