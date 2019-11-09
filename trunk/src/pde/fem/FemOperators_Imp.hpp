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
#include "FemHexahedron.hpp"
#include "FemQuadrilateral.hpp"
#include "FemTetrahedron.hpp"
#include "FemTriangle.hpp"
#include "FemTriangularPrism.hpp"

using namespace ENigMA::fem;

#include "flow/FemFlowOperators_Imp.hpp"
#include "generic/FemGenericOperators_Imp.hpp"
#include "structural/FemStructuralOperators_Imp.hpp"
#include "thermal/FemThermalOperators_Imp.hpp"

namespace ENigMA {

namespace pde {

    namespace fem {

        template <typename Real>
        void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
        {

            if (aField.simulationType() == ST_GENERIC)
                generic::ddt(aSystem, aField);
            else if (aField.simulationType() == ST_THERMAL)
                thermal::ddt(aSystem, aField);
            else if (aField.simulationType() == ST_STRUCTURAL)
                structural::ddt(aSystem, aField);
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
            else if (aField.simulationType() == ST_STRUCTURAL)
                structural::laplacian(aSystem, aField);
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
            else if (aField.simulationType() == ST_STRUCTURAL)
                structural::divergence(aSystem, aField);
        }

        template <typename Real>
        void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField1, CPdeField<Real>& aField2, Real dt)
        {

            if (aField1.simulationType() == ST_FLOW)
                flow::divergence(aSystem, aField1, aField2, dt);
        }

        template <typename Real>
        void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField1, CPdeField<Real>& aField2, CPdeField<Real>& aField3, Real dt)
        {

            if (aField1.simulationType() == ST_FLOW)
                flow::divergence(aSystem, aField1, aField2, aField3, dt);
        }

        template <typename Real>
        void gradient(CSleSystem<Real>& aSystem, CPdeField<Real>& aField, const EComponent aComponent)
        {

            if (aField.simulationType() == ST_FLOW)
                flow::gradient(aSystem, aField, aComponent);
        }

        template <typename Real>
        void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, const Real aSource)
        {

            if (aField.simulationType() == ST_GENERIC)
                generic::source(aVectorB, aField, aSource);
            else if (aField.simulationType() == ST_THERMAL)
                thermal::source(aVectorB, aField, aSource);
            else if (aField.simulationType() == ST_STRUCTURAL)
                structural::source(aVectorB, aField, aSource);
            else if (aField.simulationType() == ST_FLOW)
                flow::source(aVectorB, aField, aSource);
        }
    }
}
}
