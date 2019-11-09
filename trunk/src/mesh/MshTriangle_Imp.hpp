// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

using namespace ENigMA::geometry;

namespace ENigMA {

namespace mesh {

    template <typename Real>
    CMshTriangle<Real>::CMshTriangle()
        : m_quality(0.0)
    {
    }

    template <typename Real>
    CMshTriangle<Real>::~CMshTriangle()
    {
    }

    template <typename Real>
    void CMshTriangle<Real>::calculateQuality()
    {

        CGeoVector<Real> v0 = this->m_vertices[1] - this->m_vertices[0];
        CGeoVector<Real> v1 = this->m_vertices[2] - this->m_vertices[1];
        CGeoVector<Real> v2 = this->m_vertices[0] - this->m_vertices[2];

        Real a = v0.angle(-v2);
        Real b = v1.angle(-v0);
        Real c = v2.angle(-v1);

        Real d = sin(a) + sin(b) + sin(c);

        if (d > 0)
            m_quality = 4.0 * sin(a) * sin(b) * sin(c) / d;
        else
            m_quality = 0.0;
    }

    template <typename Real>
    Real CMshTriangle<Real>::quality()
    {

        return m_quality;
    }
}
}
