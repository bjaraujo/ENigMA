// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA {

namespace fem {

    template <typename Real>
    CFemBeam<Real, 2, 1, 1>::CFemBeam()
    {

        this->ddt.resize(2, 2);
        this->laplacian.resize(2, 2);
        this->divergence.resize(2, 2);

        this->source.resize(2);
        this->source.setZero();

        this->m_sectionArea = 1.0;
        this->m_perimeter = 1.0;
    }

    template <typename Real>
    CFemBeam<Real, 2, 1, 1>::~CFemBeam()
    {
    }

    template <typename Real>
    inline void CFemBeam<Real, 2, 1, 1>::setTransientTerm()
    {

        this->ddt << 1.0 / 3.0, 1.0 / 6.0,
            1.0 / 6.0, 1.0 / 3.0;

        this->ddt *= this->m_length * this->m_sectionArea;

        if (this->m_dt > 0.0)
            this->ddt /= this->m_dt;
    }

    template <typename Real>
    inline void CFemBeam<Real, 2, 1, 1>::setDiffusionTerm()
    {

        CFemElement<Real>::laplacian << -1, +1,
            +1, -1;

        if (this->m_length > 0)
            CFemElement<Real>::laplacian *= CFemElement<Real>::diffusionCoefficient() * this->m_sectionArea / this->m_length;
    }

    template <typename Real>
    void CFemBeam<Real, 2, 1, 1>::setConvectiveTerm()
    {

        CFemElement<Real>::divergence << -1, +1,
            -1, +1;

        CFemElement<Real>::divergence *= CFemElement<Real>::convectionCoefficient() * 0.5;
    }

    template <typename Real>
    void CFemBeam<Real, 2, 1, 1>::rebuild()
    {
    }

    template <typename Real>
    void CFemBeam<Real, 2, 1, 1>::update()
    {

        // Calculate geometrical properties
        rebuild();

        // Diffusion term
        setDiffusionTerm();

        // Convection term
        setConvectiveTerm();

        // Transient term
        if (CFemElement<Real>::transient())
            setTransientTerm();
    }

    template <typename Real>
    void CFemBeam<Real, 2, 1, 1>::setSourceOnNode(const Integer aNodeIndex, const Real aValue)
    {

        CFemElement<Real>::source(aNodeIndex) += aValue;
    }
}
}
