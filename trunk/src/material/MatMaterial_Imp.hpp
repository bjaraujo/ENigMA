// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {

namespace material {

    template <typename Real>
    CMatMaterial<Real>::CMatMaterial()
    {
    }

    template <typename Real>
    CMatMaterial<Real>::~CMatMaterial()
    {
    }

    template <typename Real>
    void CMatMaterial<Real>::addProperty(EPropertyType aPropertyType, const Real aPropertyValue)
    {

        m_property[aPropertyType] = aPropertyValue;
    }

    template <typename Real>
    void CMatMaterial<Real>::addProperty(EPropertyType aPropertyType, CAnaFunction<Real>& aPropertyFunction)
    {

        m_propertyFunc[aPropertyType] = aPropertyFunction;
    }

    template <typename Real>
    Real CMatMaterial<Real>::propertyValue(EPropertyType aPropertyType)
    {

        return m_property[aPropertyType];
    }

    template <typename Real>
    Real CMatMaterial<Real>::propertyValue(EPropertyType aPropertyType, const std::string aVariable1, Real& aValue1)
    {

        return m_propertyFunc[aPropertyType].evaluate(aVariable1, aValue1);
    }

    template <typename Real>
    Real CMatMaterial<Real>::propertyValue(EPropertyType aPropertyType, const std::string aVariable1, Real& aValue1, std::string aVariable2, Real& aValue2)
    {

        return m_propertyFunc[aPropertyType].evaluate(aVariable1, aValue1, aVariable2, aValue2);
    }
}
}
