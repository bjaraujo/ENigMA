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

namespace geometry {

    template <typename Real>
    CGeoVolume<Real>::CGeoVolume()
        : m_surfaceArea(0.0)
        , m_volume(0.0)
        , m_bSurfaceArea(false)
        , m_bVolume(false)
        , m_bBoundingBox(false)
    {
    }

    template <typename Real>
    CGeoVolume<Real>::~CGeoVolume()
    {
    }

    template <typename Real>
    Real& CGeoVolume<Real>::volume()
    {

        return m_volume;
    }

    template <typename Real>
    Real& CGeoVolume<Real>::surfaceArea()
    {

        return m_surfaceArea;
    }

    template <typename Real>
    CGeoBoundingBox<Real>& CGeoVolume<Real>::boundingBox()
    {

        return m_boundingBox;
    }
}
}
