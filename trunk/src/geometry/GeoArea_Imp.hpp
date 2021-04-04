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
    CGeoArea<Real>::CGeoArea()
        : m_area(0.0)
        , m_bNormal(false)
        , m_bArea(false)
        , m_bBoundingBox(false)
    {
    }

    template <typename Real>
    CGeoArea<Real>::~CGeoArea()
    {
    }

    template <typename Real>
    CGeoNormal<Real>& CGeoArea<Real>::normal()
    {
        return m_normal;
    }

    template <typename Real>
    Real CGeoArea<Real>::area()
    {
        return m_area;
    }

    template <typename Real>
    CGeoBoundingBox<Real>& CGeoArea<Real>::boundingBox()
    {
        return m_boundingBox;
    }
}
}
