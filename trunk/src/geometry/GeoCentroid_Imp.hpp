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
    CGeoCentroid<Real>::CGeoCentroid()
        : m_bCentroid(false)
    {
        m_centroid << 0.0, 0.0, 0.0;
    }

    template <typename Real>
    CGeoCentroid<Real>::~CGeoCentroid()
    {
    }

    template <typename Real>
    CGeoCoordinate<Real> CGeoCentroid<Real>::centroid() const
    {
        return m_centroid;
    }
}
}
