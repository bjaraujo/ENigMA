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
    CGeoLength<Real>::CGeoLength()
        : m_length(0.0)
        , m_bLength(false)
        , m_bBoundingBox(false)
    {
    }

    template <typename Real>
    CGeoLength<Real>::~CGeoLength()
    {
    }

    template <typename Real>
    Real CGeoLength<Real>::length()
    {

        return this->m_length;
    }

    template <typename Real>
    CGeoBoundingBox<Real>& CGeoLength<Real>::boundingBox()
    {

        return this->m_boundingBox;
    }
}
}
