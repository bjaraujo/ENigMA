// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "AnaFunction.hpp"

using namespace ENigMA::analytical;

namespace ENigMA {

namespace material {

    enum EPropertyType {
        PT_DENSITY = 0,
        PT_VISCOSITY,
        PT_THERMAL_CONDUCTIVITY,
        PT_SPECIFIC_HEAT,
        PT_ELASTIC_MODULUS,
        PT_POISSON_COEFFICIENT
    };

    template <typename Real>
    class CMatMaterial {
    private:
        std::map<EPropertyType, Real> m_property;
        std::map<EPropertyType, CAnaFunction<Real>> m_propertyFunc;

    public:
        CMatMaterial();
        virtual ~CMatMaterial();

        void addProperty(EPropertyType aPropertyType, const Real aPropertyValue);
        void addProperty(EPropertyType aPropertyType, CAnaFunction<Real>& aPropertyFunction);

        Real propertyValue(EPropertyType aPropertyType);
        Real propertyValue(EPropertyType aPropertyType, const std::string aVariable1, Real& aValue1);
        Real propertyValue(EPropertyType aPropertyType, const std::string aVariable1, Real& aValue1, std::string aVariable2, Real& aValue2);
    };
}
}

#include "MatMaterial_Imp.hpp"
