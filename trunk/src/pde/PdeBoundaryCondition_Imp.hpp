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

namespace pde {

    template <typename Real>
    CPdeBoundaryCondition<Real>::CPdeBoundaryCondition(EBoundaryConditionType aBoundaryConditionType) :
        m_boundaryConditionType(aBoundaryConditionType),
        m_boundaryLocation(BL_NODE)
    {
    }

    template <typename Real>
    CPdeBoundaryCondition<Real>::CPdeBoundaryCondition() : 
        m_boundaryConditionType(BT_GENERIC_FIXED_VALUE),
        m_boundaryLocation(BL_NODE)
    {
    }

    template <typename Real>
    CPdeBoundaryCondition<Real>::~CPdeBoundaryCondition()
    {
    }

    template <typename Real>
    void CPdeBoundaryCondition<Real>::setBoundaryConditionType(EBoundaryConditionType aBoundaryConditionType)
    {

        m_boundaryConditionType = aBoundaryConditionType;
    }

    template <typename Real>
    EBoundaryConditionType CPdeBoundaryCondition<Real>::boundaryConditionType()
    {

        return m_boundaryConditionType;
    }

    template <typename Real>
    void CPdeBoundaryCondition<Real>::setLocation(EBoundaryLocation aBoundaryLocation)
    {

        m_boundaryLocation = aBoundaryLocation;
    }

    template <typename Real>
    EBoundaryLocation CPdeBoundaryCondition<Real>::location()
    {

        return m_boundaryLocation;
    }

    template <typename Real>
    void CPdeBoundaryCondition<Real>::addCondition(EConditionType aConditionType, const Real aPropertyValue)
    {

        m_condition[aConditionType] = aPropertyValue;
    }

    template <typename Real>
    void CPdeBoundaryCondition<Real>::addCondition(EConditionType aConditionType, CAnaFunction<Real>& aPropertyFunction)
    {

        m_conditionFunc[aConditionType] = aPropertyFunction;
    }

    template <typename Real>
    Real CPdeBoundaryCondition<Real>::conditionValue(EConditionType aConditionType)
    {

        return m_condition[aConditionType];
    }

    template <typename Real>
    Real CPdeBoundaryCondition<Real>::conditionValue(EConditionType aConditionType, std::string aVariable1, Real& aValue1)
    {

        return m_conditionFunc[aConditionType].evaluate(aVariable1, aValue1);
    }
}
}
