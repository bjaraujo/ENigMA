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

    template <typename Real>
    CFemElement<Real>::CFemElement()
    {

        m_dt = 1.0;
        m_diffusionCoefficient = 1.0;
        m_convectionCoefficient = 1.0;

        m_transient = false;
    }

    template <typename Real>
    CFemElement<Real>::~CFemElement()
    {
    }

    template <typename Real>
    void CFemElement<Real>::setDt(const Real aValue)
    {

        m_dt = aValue;
    }

    template <typename Real>
    Real CFemElement<Real>::dt() const
    {

        return m_dt;
    }

    template <typename Real>
    void CFemElement<Real>::setDiffusionCoefficient(const Real aValue)
    {

        m_diffusionCoefficient = aValue;
    }

    template <typename Real>
    Real CFemElement<Real>::diffusionCoefficient() const
    {

        return m_diffusionCoefficient;
    }

    template <typename Real>
    void CFemElement<Real>::setConvectionCoefficient(const Real aValue)
    {

        m_convectionCoefficient = aValue;
    }

    template <typename Real>
    Real CFemElement<Real>::convectionCoefficient() const
    {

        return m_convectionCoefficient;
    }

    template <typename Real>
    void CFemElement<Real>::setTransient(const bool aValue)
    {

        m_transient = aValue;
    }

    template <typename Real>
    bool CFemElement<Real>::transient() const
    {

        return m_transient;
    }
}
}
