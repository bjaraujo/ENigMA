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
namespace fem {
    namespace structural {
        template <typename Real>
        CFemStructuralElement<Real>::CFemStructuralElement()
            : m_coeffPoisson(1.0)
            , m_elasticModulus(1.0)
        {
        }

        template <typename Real>
        CFemStructuralElement<Real>::~CFemStructuralElement()
        {
        }

        template <typename Real>
        void CFemStructuralElement<Real>::setCoeffPoisson(const Real aValue)
        {
            m_coeffPoisson = aValue;
        }

        template <typename Real>
        Real CFemStructuralElement<Real>::coeffPoisson() const
        {
            return m_coeffPoisson;
        }

        template <typename Real>
        void CFemStructuralElement<Real>::setElasticModulus(const Real aValue)
        {
            m_elasticModulus = aValue;
        }

        template <typename Real>
        Real CFemStructuralElement<Real>::elasticModulus() const
        {
            return m_elasticModulus;
        }
    }
}
}
