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

#include "../../src/fem/FemTetrahedron.hpp"
#include "FemThermalElement.hpp"

namespace ENigMA
{

    namespace fem
    {

        namespace thermal
        {

            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemLinearTemperatureTetrahedron : public CFemThermalElement<Real>, public CFemTetrahedron<Real, NbNodes, Dof, Order>
            {
            };

            template <typename Real>
            class CFemLinearTemperatureTetrahedron<Real, 4, 1, 1> : public CFemThermalElement<Real>, public CFemTetrahedron<Real, 4, 1, 1>
            {
            public:
                
                CFemLinearTemperatureTetrahedron();
                ~CFemLinearTemperatureTetrahedron();

                void setConvectionOnFace(const Integer aFaceIndex, const Real h, const Real Tinf);
                void setConvectionOnFace(const Integer aFaceIndex, const Real e, const Real teta, const Real Tinf);

            };

        }

    }

}

#include "FemLinearTemperatureTetrahedron_Imp.hpp"


