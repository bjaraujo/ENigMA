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
    CMshQuadrilateral<Real>::CMshQuadrilateral()
        : m_quality(0.0)
    {
    }

    template <typename Real>
    CMshQuadrilateral<Real>::~CMshQuadrilateral()
    {
    }

    template <typename Real>
    void CMshQuadrilateral<Real>::calculateQuality()
    {

        this->calculateArea();

        Real Ai = CGeoQuadrilateral<Real>::area();

        Real ai = (this->m_vertices[1] - this->m_vertices[0]).squaredNorm();
        Real bi = (this->m_vertices[2] - this->m_vertices[1]).squaredNorm();
        Real ci = (this->m_vertices[3] - this->m_vertices[2]).squaredNorm();
        Real di = (this->m_vertices[0] - this->m_vertices[3]).squaredNorm();

        m_quality = 4.0 * Ai / (ai + bi + ci + di);
    }

    template <typename Real>
    Real CMshQuadrilateral<Real>::quality()
    {

        return m_quality;
    }
}
}
