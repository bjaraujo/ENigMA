// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

namespace ENigMA {

namespace mesh {

    template <typename Real>
    CMshTetrahedron<Real>::CMshTetrahedron()
        : m_quality(0.0)
    {
    }

    template <typename Real>
    CMshTetrahedron<Real>::~CMshTetrahedron()
    {
    }

    template <typename Real>
    void CMshTetrahedron<Real>::calculateQuality()
    {

        this->calculateSurfaceArea();
        this->calculateVolume();

        Real S = this->m_surfaceArea;
        Real V = this->m_volume;

        Real l1 = (this->m_vertices[1] - this->m_vertices[0]).norm();
        Real l2 = (this->m_vertices[2] - this->m_vertices[1]).norm();
        Real l3 = (this->m_vertices[0] - this->m_vertices[2]).norm();
        Real l4 = (this->m_vertices[3] - this->m_vertices[0]).norm();
        Real l5 = (this->m_vertices[3] - this->m_vertices[1]).norm();
        Real l6 = (this->m_vertices[3] - this->m_vertices[2]).norm();

        Real maxL = std::max(l1, std::max(l2, std::max(l3, std::max(l4, std::max(l5, l6)))));

        if (fabs(S * maxL) > 0.0)
            m_quality = 6.0 * sqrt(6.0) * fabs(V) / fabs(S * maxL);
        else
            m_quality = 0.0;
    }

    template <typename Real>
    Real CMshTetrahedron<Real>::quality()
    {

        return m_quality;
    }
}
}
