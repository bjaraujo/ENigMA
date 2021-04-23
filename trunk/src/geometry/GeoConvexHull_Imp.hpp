// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        CGeoConvexHull<Real>::CGeoConvexHull()
        {
        }

        template <typename Real>
        CGeoConvexHull<Real>::CGeoConvexHull(const std::vector<CGeoCoordinate<Real>>& sVertices)
        {
            this->set(sVertices);
        }

        template <typename Real>
        void CGeoConvexHull<Real>::set(const std::vector<CGeoCoordinate<Real>>& sVertices)
        {
            if (sVertices.size() < 3)
                return;

            CGeoVector<Real> v0 = sVertices.at(1) - sVertices.at(0);
            CGeoVector<Real> v1 = sVertices.at(2) - sVertices.at(0);

            CGeoVector<Real> vx = v0;
            vx.normalize();

            CGeoVector<Real> vz = vx.cross(v1);
            vz.normalize();

            CGeoVector<Real> vy = vz.cross(vx);
            vy.normalize();

            CGeoCoordinateSystem<Real> aCoordinateSystem;
            aCoordinateSystem.col(0) << vx;
            aCoordinateSystem.col(1) << vy;
            aCoordinateSystem.col(2) << vz;

            // Project to 2D plane
            std::vector<SVertex<Real>> sExtVertices;

            for (Integer i = 0; i < static_cast<Integer>(sVertices.size()); ++i)
            {
                SVertex<Real> aVertex;

                aVertex.v2D = sVertices.at(i);
                aVertex.v3D = sVertices.at(i);

                aVertex.v2D.transform(aCoordinateSystem);

                sExtVertices.push_back(aVertex);
            }

            Integer n = static_cast<Integer>(sExtVertices.size());
            Integer k = 0;

            std::vector<SVertex<Real>> sNewExtVertices(2 * n);

            // Sort points lexicographically
            std::sort(sExtVertices.begin(), sExtVertices.end());

            // Build lower hull
            for (Integer i = 0; i < n; ++i)
            {
                while (k >= 2)
                {
                    Real c = (sNewExtVertices[k - 2].v2D - sNewExtVertices[k - 1].v2D).cross(sNewExtVertices[k - 2].v2D - sExtVertices[i].v2D).z();

                    if (c <= 0)
                        break;

                    k--;
                }

                sNewExtVertices[k++] = sExtVertices[i];
            }

            // Build upper hull
            for (Integer i = n - 2, t = k + 1; i >= 0; i--)
            {
                while (k >= t)
                {
                    Real c = (sNewExtVertices[k - 2].v2D - sNewExtVertices[k - 1].v2D).cross(sNewExtVertices[k - 2].v2D - sExtVertices[i].v2D).z();

                    if (c <= 0)
                        break;

                    k--;
                }

                sNewExtVertices[k++] = sExtVertices[i];
            }

            sNewExtVertices.resize(k - 1);

            for (Integer i = 0; i < static_cast<Integer>(sNewExtVertices.size()); ++i)
                m_vertices.push_back(sNewExtVertices[i].v3D);
        }

        template <typename Real>
        CGeoConvexHull<Real>::~CGeoConvexHull()
        {
            reset();
        }

        template <typename Real>
        void CGeoConvexHull<Real>::reset()
        {
            m_vertices.clear();
        }

        template <typename Real>
        Integer CGeoConvexHull<Real>::nbVertices() const
        {
            return static_cast<Integer>(m_vertices.size());
        }

        template <typename Real>
        void CGeoConvexHull<Real>::addVertex(CGeoCoordinate<Real>& aVertex)
        {
            m_vertices.push_back(aVertex);
        }

        template <typename Real>
        CGeoCoordinate<Real>& CGeoConvexHull<Real>::vertex(const Integer aVertexIndex)
        {
            return m_vertices[aVertexIndex % m_vertices.size()];
        }
    }
}
