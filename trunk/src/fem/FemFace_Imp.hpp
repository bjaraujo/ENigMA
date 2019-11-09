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
    CFemFace<Real>::CFemFace()
    {
    }

    template <typename Real>
    CFemFace<Real>::~CFemFace()
    {
    }

    template <typename Real>
    void CFemFace<Real>::setThickness(const Real aValue)
    {

        m_thickness = aValue;
    }

    template <typename Real>
    Real CFemFace<Real>::thickness() const
    {

        return m_thickness;
    }
}
}
