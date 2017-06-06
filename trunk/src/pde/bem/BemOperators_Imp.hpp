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

#include "../../src/bem/BemTriangle.hpp"

using namespace ENigMA::bem;

#include "thermal/BemThermalOperators_Imp.hpp"

namespace ENigMA
{

    namespace pde
    {

        namespace bem
        {

            template <typename Real>
            void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
            {

            }

            template <typename Real>
            void laplacian(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
            {

                if (aField.simulationType() == ST_THERMAL)
                    thermal::laplacian(aSystem, aField);

            }

            template <typename Real>
            void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
            {

                // TODO:

            }

            template <typename Real>
            void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, Real aSource)
            {

            }

        }

    }

}

