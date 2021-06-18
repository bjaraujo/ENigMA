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
        CGeoPolygon<Real>::CGeoPolygon()
        {
        }

        template <typename Real>
        CGeoPolygon<Real>::~CGeoPolygon()
        {
            reset();
        }

        template <typename Real>
        void CGeoPolygon<Real>::reset()
        {
            m_polyline.reset();
        }

        template <typename Real>
        CGeoPolyline<Real>& CGeoPolygon<Real>::polyline()
        {
            return m_polyline;
        }

        template <typename Real>
        void CGeoPolygon<Real>::setPolyline(const CGeoPolyline<Real>& aPolyline)
        {
            m_polyline = aPolyline;
        }

        template <typename Real>
        void CGeoPolygon<Real>::calculateCentroid(bool bReCalculate)
        {
            if (!this->m_bCentroid || bReCalculate)
            {
                this->m_centroid << 0.0, 0.0, 0.0;

                for (Integer i = 1; i < m_polyline.nbVertices(); ++i)
                {
                    this->m_centroid += m_polyline.vertex(i);
                }

                if (m_polyline.nbVertices() > 1)
                    this->m_centroid /= (static_cast<Real>(m_polyline.nbVertices()) - 1);

                this->m_bCentroid = true;
            }
        }

        template <typename Real>
        void CGeoPolygon<Real>::calculateNormal(bool bReCalculate)
        {
            if (!this->m_bNormal || bReCalculate)
            {
                CGeoVector<Real> v0, v1;
                CGeoNormal<Real> aNormal;

                for (Integer i = 1; i < m_polyline.nbVertices() - 1; ++i)
                {
                    // Calculate normal vector of the polygon's plane
                    v0 = m_polyline.vertex(i + 0) - m_polyline.vertex(i - 1);
                    v1 = m_polyline.vertex(i + 1) - m_polyline.vertex(i - 1);

                    aNormal = v0.cross(v1);

                    if (aNormal.norm() > std::numeric_limits<Real>::epsilon())
                    {
                        aNormal.normalize();
                        break;
                    }
                }

                this->m_normal = aNormal;

                this->m_bNormal = true;
            }
        }

        template <typename Real>
        void CGeoPolygon<Real>::calculateArea(bool bReCalculate)
        {
            if (!this->m_bArea || bReCalculate)
            {
                this->m_area = 0.0;

                if (m_polyline.nbVertices() > 1)
                {
                    this->calculateNormal(bReCalculate);

                    if (this->m_normal.norm() > std::numeric_limits<Real>::epsilon())
                    {
                        Integer coord = 3;
                        if (fabs(this->m_normal.x()) > fabs(this->m_normal.y()))
                        {
                            if (fabs(this->m_normal.x()) > fabs(this->m_normal.z()))
                                coord = 1; // ignore x-coord
                            else
                                coord = 3; // ignore z-coord
                        }
                        else if (fabs(this->m_normal.y()) > fabs(this->m_normal.z()))
                            coord = 2; // ignore y-coord

                        Real _parea = 0.0;

                        for (Integer i = 1; i <= m_polyline.nbVertices(); ++i)
                        {
                            switch (coord)
                            {
                            case 1:
                                _parea += (m_polyline.vertex(i).y() * (m_polyline.vertex(i + 1).z() - m_polyline.vertex(i - 1).z()));
                                break;
                            case 2:
                                _parea += (m_polyline.vertex(i).z() * (m_polyline.vertex(i + 1).x() - m_polyline.vertex(i - 1).x()));
                                break;
                            case 3:
                                _parea += (m_polyline.vertex(i).x() * (m_polyline.vertex(i + 1).y() - m_polyline.vertex(i - 1).y()));
                                break;
                            }
                        }

                        switch (coord)
                        {
                        case 1:
                            _parea *= (this->m_normal.norm() / (2 * this->m_normal.x()));
                            break;
                        case 2:
                            _parea *= (this->m_normal.norm() / (2 * this->m_normal.y()));
                            break;
                        case 3:
                            _parea *= (this->m_normal.norm() / (2 * this->m_normal.z()));
                        }

                        this->m_area += _parea;
                    }
                }

                this->m_bArea = true;
            }
        }

        template <typename Real>
        void CGeoPolygon<Real>::calculateBoundingBox(bool bReCalculate)
        {
            // TODO:
        }

        template <typename Real>
        CGeoPolygon<Real> CGeoPolygon<Real>::clip(CGeoPlane<Real>& aPlane, const Real aTolerance)
        {
            CGeoPolygon<Real> aPolygon;

            CGeoLineList<Real> aLineList = m_polyline.clip(aPlane);

            CGeoPolyline<Real> aPolyline(aLineList, false, aTolerance);

            aPolyline.close(aTolerance);

            aPolygon.setPolyline(aPolyline);

            return aPolygon;
        }

        template <typename Real>
        void CGeoPolygon<Real>::invert()
        {
            m_polyline.invert();

            this->m_normal *= -1;
        }

        template <typename Real>
        void CGeoPolygon<Real>::close()
        {
            m_polyline.close();
        }

        template <typename Real>
        bool CGeoPolygon<Real>::isClosed()
        {
            return m_polyline.isClosed();
        }

        template <typename Real>
        bool CGeoPolygon<Real>::intersects(CGeoPlane<Real>& aPlane)
        {
            return m_polyline.intersects(aPlane);
        }

        template <typename Real>
        std::vector<CGeoTriangle<Real>> CGeoPolygon<Real>::triangulate()
        {
            std::vector<CGeoTriangle<Real>> sTriangles;
            std::vector<Integer> chain;

            for (Integer i = 1; i < m_polyline.nbVertices(); ++i)
                chain.push_back(i - 1);

            Integer p1 = 0;
            Integer p2 = 1;
            Integer p3 = 2;

            while (chain.size() > 2)
            {
                Real minAngle = +2 * 3.142;

                Integer wi = 0;

                for (Integer i = 0; i < static_cast<Integer>(chain.size()); ++i)
                {
                    Integer im, ip, ii;

                    if (i == 0)
                    {
                        ii = chain[i];
                        im = chain[chain.size() - 1];
                        ip = chain[1];
                    }
                    else
                    {
                        ii = chain[i];
                        im = chain[i - 1];
                        ip = chain[(i + 1) % chain.size()];
                    }

                    CGeoVector<Real> v1 = m_polyline.vertex(im) - m_polyline.vertex(ii);
                    CGeoVector<Real> v2 = m_polyline.vertex(ip) - m_polyline.vertex(ii);

                    Real angle = fabs(v1.angle(v2));

                    if (angle < minAngle)
                    {
                        minAngle = angle;

                        p1 = ii;
                        p2 = ip;
                        p3 = im;

                        wi = i;
                    }
                }

                CGeoTriangle<Real> aTriangle;

                aTriangle.addVertex(m_polyline.vertex(p1));
                aTriangle.addVertex(m_polyline.vertex(p2));
                aTriangle.addVertex(m_polyline.vertex(p3));

                sTriangles.push_back(aTriangle);

                chain.erase(chain.begin() + wi);
            }

            return sTriangles;
        }
    }
}
