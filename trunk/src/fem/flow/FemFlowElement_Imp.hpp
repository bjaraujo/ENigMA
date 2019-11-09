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

    namespace flow {

        template <typename Real>
        CFemFlowElement<Real>::CFemFlowElement() : m_density(1.0), m_viscosity(1.0)
        {
        }

        template <typename Real>
        CFemFlowElement<Real>::~CFemFlowElement()
        {
        }

        template <typename Real>
        void CFemFlowElement<Real>::setDensity(const Real aValue)
        {

            m_density = aValue;
        }

        template <typename Real>
        Real CFemFlowElement<Real>::density() const
        {

            return m_density;
        }

        template <typename Real>
        void CFemFlowElement<Real>::setViscosity(const Real aValue)
        {

            m_viscosity = aValue;
        }

        template <typename Real>
        Real CFemFlowElement<Real>::viscosity() const
        {

            return m_viscosity;
        }
    }
}
}
