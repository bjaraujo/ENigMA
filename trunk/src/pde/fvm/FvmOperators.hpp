// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "PdeField.hpp"
#include "SleSystem.hpp"

using namespace ENigMA::pde;
using namespace ENigMA::sle;

namespace ENigMA
{

    namespace pde
    {

        namespace fvm
        {

            template <typename Real>
            void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField);

            template <typename Real>
            void laplacian(CSleSystem<Real>& aSystem, CPdeField<Real>& aField);

            template <typename Real>
            void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField);

            template <typename Real>
            void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, Real aSource);

        }

    }

}

#include "FvmOperators_Imp.hpp"

