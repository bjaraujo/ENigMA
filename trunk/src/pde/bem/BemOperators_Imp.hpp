// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "BemTriangle.hpp"
#include "BemThermalOperators_Imp.hpp"

using namespace ENigMA::bem;

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

