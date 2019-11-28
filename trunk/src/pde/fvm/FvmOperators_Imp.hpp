// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FvmControlVolume.hpp"

#include "flow/FvmFlowOperators_Imp.hpp"
#include "generic/FvmGenericOperators_Imp.hpp"
#include "thermal/FvmThermalOperators_Imp.hpp"

using namespace ENigMA::fvm;

namespace ENigMA {
namespace pde {
    namespace fvm {
        template <typename Real>
        void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
        {
            if (aField.simulationType() == ST_GENERIC)
                generic::ddt(aSystem, aField);
            else if (aField.simulationType() == ST_THERMAL)
                thermal::ddt(aSystem, aField);
            else if (aField.simulationType() == ST_FLOW)
                flow::ddt(aSystem, aField);
        }

        template <typename Real>
        void laplacian(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
        {
            if (aField.simulationType() == ST_GENERIC)
                generic::laplacian(aSystem, aField);
            else if (aField.simulationType() == ST_THERMAL)
                thermal::laplacian(aSystem, aField);
            else if (aField.simulationType() == ST_FLOW)
                flow::laplacian(aSystem, aField);
        }

        template <typename Real>
        void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
        {
            if (aField.simulationType() == ST_GENERIC)
                generic::divergence(aSystem, aField);
            else if (aField.simulationType() == ST_THERMAL)
                thermal::divergence(aSystem, aField);
            else if (aField.simulationType() == ST_FLOW)
                flow::divergence(aSystem, aField);
        }

        template <typename Real>
        void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, Real aSource)
        {
            if (aField.simulationType() == ST_GENERIC)
                generic::source(aVectorB, aField, aSource);
            else if (aField.simulationType() == ST_THERMAL)
                thermal::source(aVectorB, aField, aSource);
            else if (aField.simulationType() == ST_FLOW)
                flow::source(aVectorB, aField, aSource);
        }
    }
}
}
