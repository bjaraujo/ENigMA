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

#include "../../src/fem/FemTriangle.hpp"
#include "FemThermalElement.hpp"

namespace ENigMA
{

    namespace fem
    {

        namespace thermal
        {

            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemLinearTemperatureTriangle : public CFemThermalElement<Real>, public CFemTriangle<Real, NbNodes, Dof, Order>
            {
            };

            template <typename Real>
            class CFemLinearTemperatureTriangle<Real, 3, 1, 1> : public CFemThermalElement<Real>, public CFemTriangle<Real, 3, 1, 1>
            {
            public:
                
                CFemLinearTemperatureTriangle();
                ~CFemLinearTemperatureTriangle();

                void setDiffusionTerm();

                void setConvectionOnEdge(const Integer anEdgeIndex, const Real h, const Real Tinf);
                void setConvectionOnEdge(const Integer anEdgeIndex, const Real e, const Real teta, const Real Tinf);

            };

        }

    }

}

#include "FemLinearTemperatureTriangle_Imp.hpp"


