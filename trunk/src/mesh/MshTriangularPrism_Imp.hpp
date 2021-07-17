// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        CMshTriangularPrism<Real>::CMshTriangularPrism()
            : m_quality(0.0)
        {
        }

        template <typename Real>
        CMshTriangularPrism<Real>::~CMshTriangularPrism()
        {
        }

        template <typename Real>
        void CMshTriangularPrism<Real>::calculateQuality()
        {
            m_quality = 0.0;
            
            {
                CGeoVector<Real> v0 = this->m_vertices[1] - this->m_vertices[0];
                CGeoVector<Real> v1 = this->m_vertices[2] - this->m_vertices[1];
                CGeoVector<Real> v2 = this->m_vertices[0] - this->m_vertices[2];

                Real a = v0.angle(-v2);
                Real b = v1.angle(-v0);
                Real c = v2.angle(-v1);

                Real d = sin(a) + sin(b) + sin(c);

                if (d > std::numeric_limits<Real>::epsilon())
                    m_quality += 2.0 * sin(a) * sin(b) * sin(c) / d;                              
            }

            {
                CGeoVector<Real> v0 = this->m_vertices[4] - this->m_vertices[3];
                CGeoVector<Real> v1 = this->m_vertices[5] - this->m_vertices[4];
                CGeoVector<Real> v2 = this->m_vertices[3] - this->m_vertices[5];

                Real a = v0.angle(-v2);
                Real b = v1.angle(-v0);
                Real c = v2.angle(-v1);

                Real d = sin(a) + sin(b) + sin(c);

                if (d > std::numeric_limits<Real>::epsilon())
                    m_quality += 2.0 * sin(a) * sin(b) * sin(c) / d;                              
            }
        }

        template <typename Real>
        Real CMshTriangularPrism<Real>::quality() const
        {
            return m_quality;
        }
    }
}
