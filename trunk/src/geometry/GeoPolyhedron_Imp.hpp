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
    CGeoPolyhedron<Real>::CGeoPolyhedron()
    {
    }

    template <typename Real>
    CGeoPolyhedron<Real>::CGeoPolyhedron(CGeoTetrahedron<Real>& aTetrahedron)
    {
        CGeoPolyline<Real> aPolyline;
        CGeoPolygon<Real> aPolygon;

        // face 1
        aPolyline.addVertex(aTetrahedron.vertex(0));
        aPolyline.addVertex(aTetrahedron.vertex(2));
        aPolyline.addVertex(aTetrahedron.vertex(1));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(0, aPolygon);
        aPolyline.reset();
        aPolygon.reset();

        // face 2
        aPolyline.addVertex(aTetrahedron.vertex(0));
        aPolyline.addVertex(aTetrahedron.vertex(1));
        aPolyline.addVertex(aTetrahedron.vertex(3));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(1, aPolygon);
        aPolyline.reset();
        aPolygon.reset();

        // face 3
        aPolyline.addVertex(aTetrahedron.vertex(0));
        aPolyline.addVertex(aTetrahedron.vertex(3));
        aPolyline.addVertex(aTetrahedron.vertex(2));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(2, aPolygon);
        aPolyline.reset();
        aPolygon.reset();

        // face 4
        aPolyline.addVertex(aTetrahedron.vertex(1));
        aPolyline.addVertex(aTetrahedron.vertex(2));
        aPolyline.addVertex(aTetrahedron.vertex(3));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(3, aPolygon);
        aPolyline.reset();
        aPolygon.reset();
    }

    template <typename Real>
    CGeoPolyhedron<Real>::CGeoPolyhedron(CGeoTriangularPrism<Real>& aTriangularPrism)
    {
        // TODO:
    }

    template <typename Real>
    CGeoPolyhedron<Real>::CGeoPolyhedron(CGeoHexahedron<Real>& aHexahedron)
    {
        CGeoPolyline<Real> aPolyline;
        CGeoPolygon<Real> aPolygon;

        // face 1
        aPolyline.addVertex(aHexahedron.vertex(0));
        aPolyline.addVertex(aHexahedron.vertex(3));
        aPolyline.addVertex(aHexahedron.vertex(2));
        aPolyline.addVertex(aHexahedron.vertex(1));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(0, aPolygon);
        aPolyline.reset();
        aPolygon.reset();

        // face 2
        aPolyline.addVertex(aHexahedron.vertex(4));
        aPolyline.addVertex(aHexahedron.vertex(5));
        aPolyline.addVertex(aHexahedron.vertex(6));
        aPolyline.addVertex(aHexahedron.vertex(7));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(1, aPolygon);
        aPolyline.reset();
        aPolygon.reset();

        // face 3
        aPolyline.addVertex(aHexahedron.vertex(4));
        aPolyline.addVertex(aHexahedron.vertex(0));
        aPolyline.addVertex(aHexahedron.vertex(1));
        aPolyline.addVertex(aHexahedron.vertex(5));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(2, aPolygon);
        aPolyline.reset();
        aPolygon.reset();

        // face 4
        aPolyline.addVertex(aHexahedron.vertex(7));
        aPolyline.addVertex(aHexahedron.vertex(6));
        aPolyline.addVertex(aHexahedron.vertex(2));
        aPolyline.addVertex(aHexahedron.vertex(3));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(3, aPolygon);
        aPolyline.reset();
        aPolygon.reset();

        // face 5
        aPolyline.addVertex(aHexahedron.vertex(5));
        aPolyline.addVertex(aHexahedron.vertex(1));
        aPolyline.addVertex(aHexahedron.vertex(2));
        aPolyline.addVertex(aHexahedron.vertex(6));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(4, aPolygon);
        aPolyline.reset();
        aPolygon.reset();

        // face 6
        aPolyline.addVertex(aHexahedron.vertex(0));
        aPolyline.addVertex(aHexahedron.vertex(4));
        aPolyline.addVertex(aHexahedron.vertex(7));
        aPolyline.addVertex(aHexahedron.vertex(3));
        aPolyline.close();

        aPolygon.setPolyline(aPolyline);
        this->addPolygon(5, aPolygon);
        aPolyline.reset();
        aPolygon.reset();
    }

    template <typename Real>
    CGeoPolyhedron<Real>::~CGeoPolyhedron()
    {
        reset();
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::reset()
    {
        m_polygons.clear();
        m_polygonIds.clear();
    }

    template <typename Real>
    Integer CGeoPolyhedron<Real>::polygonId(const Integer aPolygonIndex)
    {
        return m_polygonIds[aPolygonIndex];
    }

    template <typename Real>
    CGeoPolygon<Real>& CGeoPolyhedron<Real>::polygon(const Integer aPolygonId)
    {
        return m_polygons[aPolygonId];
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::addPolygon(const Integer aPolygonId, CGeoPolygon<Real>& aPolygon)
    {
        m_polygons[aPolygonId] = aPolygon;
        m_polygonIds.push_back(aPolygonId);
    }

    template <typename Real>
    Integer CGeoPolyhedron<Real>::nbPolygons() const
    {
        return static_cast<Integer>(m_polygonIds.size());
    }

    template <typename Real>
    bool CGeoPolyhedron<Real>::containsPolygon(const Integer aPolygonId)
    {
        return std::find(m_polygonIds.begin(), m_polygonIds.end(), aPolygonId) != m_polygonIds.end();
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::triangulate()
    {
        // TODO
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::close(CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoNormal<Real>& aNormal, const Real aTolerance)
    {
        // Discover lines
        CGeoLineList<Real> aLineList;

        for (Integer i = 0; i < static_cast<Integer>(m_polygonIds.size()); ++i) {
            CGeoPolyline<Real> aPolyline = m_polygons[m_polygonIds[i]].polyline();

            for (Integer j = 0; j < aPolyline.nbLines(); ++j) {
                CGeoLine<Real> aLine = aPolyline.line(j);

                aLine.calculateLength();

                if (aLine.length() > aTolerance)
                    aLineList.addLine(aLine);
            }
        }

        // Edges
        CGeoLineList<Real> anEdgeList;

        for (Integer i = 0; i < aLineList.nbLines(); ++i) {
            bool pair = false;

            for (Integer j = 0; j < aLineList.nbLines(); ++j) {
                if (i != j) {
                    Real d1, d2;

                    // Normal
                    d1 = (aLineList.line(i).startPoint() - aLineList.line(j).startPoint()).norm();
                    d2 = (aLineList.line(i).endPoint() - aLineList.line(j).endPoint()).norm();

                    if (d1 <= aTolerance && d2 <= aTolerance) {
                        pair = true;
                        break;
                    }

                    // Inverted
                    d1 = (aLineList.line(i).startPoint() - aLineList.line(j).endPoint()).norm();
                    d2 = (aLineList.line(i).endPoint() - aLineList.line(j).startPoint()).norm();

                    if (d1 <= aTolerance && d2 <= aTolerance) {
                        pair = true;
                        break;
                    }
                }
            }

            if (!pair) {
                anEdgeList.addLine(aLineList.line(i));
            }
        }

        // Add new polygon
        if (anEdgeList.nbLines() > 2) {
            CGeoPolyline<Real> aPolyline(anEdgeList, true, aTolerance);

            aNewPolygon.setPolyline(aPolyline);

            aNewPolygon.calculateNormal(true);

            if (aNormal.dot(aNewPolygon.normal()) < 0)
                aNewPolygon.invert();

            this->addPolygon(aNewPolygonId, aNewPolygon);
        }
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::calculateCentroid(bool bReCalculate)
    {
        if (!m_bCentroid || bReCalculate) {
            m_centroid << 0.0, 0.0, 0.0;

            for (Integer i = 0; i < static_cast<Integer>(m_polygonIds.size()); ++i) {
                m_polygons[m_polygonIds[i]].calculateCentroid();
                m_centroid += m_polygons[m_polygonIds[i]].centroid();
            }

            if (m_polygons.size() > 0)
                m_centroid /= static_cast<Real>(m_polygonIds.size());

            m_bCentroid = true;
        }
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::calculateSurfaceArea(bool bReCalculate)
    {
        if (!m_bSurfaceArea || bReCalculate) {
            m_surfaceArea = 0.0;

            for (Integer i = 0; i < static_cast<Integer>(m_polygonIds.size()); ++i) {
                m_polygons[m_polygonIds[i]].calculateArea();
                m_surfaceArea += m_polygons[m_polygonIds[i]].area();
            }

            m_bSurfaceArea = true;
        }
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::calculateVolume(bool bReCalculate)
    {
        if (!m_bVolume || bReCalculate) {
            // see: http://en.wikipedia.org/wiki/Polyhedron#Volume
            // see: http://g3d.sourceforge.net/

            m_volume = 0.0;

            if (m_polygons.size() >= 4) {
                if (m_polygons[m_polygonIds[0]].polyline().nbVertices() > 0) {
                    CGeoCoordinate<Real> v0 = m_polygons[m_polygonIds[0]].polyline().vertex(0);

                    for (Integer i = 1; i < static_cast<Integer>(m_polygonIds.size()); ++i) {
                        m_polygons[m_polygonIds[i]].calculateArea();

                        if (m_polygons[m_polygonIds[i]].polyline().nbVertices() > 0)
                            m_volume += (m_polygons[m_polygonIds[i]].polyline().vertex(0) - v0).dot(m_polygons[m_polygonIds[i]].normal()) * m_polygons[m_polygonIds[i]].area();
                    }
                }
            }

            m_volume /= 3.0;
            m_bVolume = true;
        }
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::calculateBoundingBox(bool bReCalculate)
    {
        if (!m_bBoundingBox || bReCalculate) {
            CGeoVolume<Real>::boundingBox().reset();

            for (Integer i = 0; i < static_cast<Integer>(m_polygonIds.size()); ++i) {
                for (Integer j = 0; j < m_polygons[m_polygonIds[i]].polyline().nbVertices(); ++j) {
                    CGeoVolume<Real>::boundingBox().addCoordinate(m_polygons[m_polygonIds[i]].polyline().vertex(j));
                }
            }

            m_bBoundingBox = true;
        }
    }

    template <typename Real>
    CGeoPolyhedron<Real> CGeoPolyhedron<Real>::clip(CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoPlane<Real>& aPlane, const Real aTolerance)
    {
        CGeoPolyhedron<Real> aPolyhedron;

        for (Integer i = 0; i < static_cast<Integer>(m_polygonIds.size()); ++i) {
            CGeoPolygon<Real> aPolygon = m_polygons[m_polygonIds[i]].clip(aPlane, aTolerance);

            if (aPolygon.polyline().nbVertices() > 3)
                aPolyhedron.addPolygon(m_polygonIds[i], aPolygon);
        }

        aPolyhedron.close(aNewPolygon, aNewPolygonId, aPlane.normal(), aTolerance);

        return aPolyhedron;
    }

    template <typename Real>
    CGeoPolyhedron<Real> CGeoPolyhedron<Real>::clip(CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoNormal<Real>& aNormal, Real& d, Real volumeFractionReq, Real& volumeFractionAct, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance)
    {
        CGeoPolyhedron<Real> aPolyhedron;

        if (volumeFractionReq <= 0.0) {
            nIterations = 0;

            volumeFractionAct = 0.0;

            return aPolyhedron;
        }

        if (volumeFractionReq >= 1.0) {
            nIterations = 0;

            // Copy this polyhedron
            for (Integer i = 0; i < static_cast<Integer>(m_polygonIds.size()); ++i) {
                aPolyhedron.addPolygon(m_polygonIds[i], m_polygons[m_polygonIds[i]]);
            }

            volumeFractionAct = 1.0;

            return aPolyhedron;
        }

        CGeoPlane<Real> aPlane(aNormal, 0.0);

        // Find interval a, b

        this->calculateVolume();
        Real vt = CGeoVolume<Real>::volume();

        Real a = +std::numeric_limits<Real>::max();
        Real b = -std::numeric_limits<Real>::max();

        for (Integer i = 0; i < static_cast<Integer>(m_polygonIds.size()); ++i) {
            for (Integer j = 0; j < m_polygons[m_polygonIds[i]].polyline().nbVertices(); ++j) {
                d = m_polygons[m_polygonIds[i]].polyline().vertex(j).dot(aNormal);

                a = std::min(a, d);
                b = std::max(b, d);
            }
        }

        // Solve using brent's method

        Real e = a;
        Real f = std::numeric_limits<Real>::max();

        Real fa = 0.0 - volumeFractionReq;
        Real fb = 1.0 - volumeFractionReq;

        Real fc = fa;
        Real s = 0.0;

        Real vs = 0.0;

        bool mflag = true;

        nIterations = nMaxIterations;

        for (Integer i = 0; i < nMaxIterations; ++i) {
            if ((fb == 0) || (fabs(fa - fb) <= aTolerance)) {
                nIterations = i + 1;
                break;
            }

            if ((fa != fc) && (fb != fc))
                // Inverse quadratic interpolation
                s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + e * fa * fb / (fc - fa) / (fc - fb);
            else
                // Secant Rule
                s = b - fb * (b - a) / (fb - fa);

            Real tmp2 = (3 * a + b) / 4;

            if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) || (mflag && (fabs(s - b) >= (fabs(b - e) / 2))) || (!mflag && (fabs(s - b) >= (fabs(e - f) / 2)))) {
                s = (a + b) * 0.5;
                mflag = true;
            } else {
                if ((mflag && (fabs(b - e) < aTolerance)) || (!mflag && (fabs(e - f) < aTolerance))) {
                    s = (a + b) / 2;
                    mflag = true;
                } else
                    mflag = false;
            }

            aPlane.setD(s);
            aPolyhedron = this->clip(aNewPolygon, aNewPolygonId, aPlane, aTolerance);

            aPolyhedron.calculateVolume(true);
            vs = aPolyhedron.volume();
            Real fs = vs / vt - volumeFractionReq;

            f = e;
            e = b;
            fc = fb;

            if (fa * fs < 0) {
                b = s;
                fb = fs;
            } else {
                a = s;
                fa = fs;
            }

            // if |f(a)| < |f(b)| then swap (a,b) end if
            if (fabs(fa) < fabs(fb)) {
                Real tmp = a;
                a = b;
                b = tmp;
                tmp = fa;
                fa = fb;
                fb = tmp;
            }
        }

        volumeFractionAct = vs / vt;

        d = s;

        return aPolyhedron;
    }

    template <typename Real>
    CGeoPolyhedron<Real> CGeoPolyhedron<Real>::cut(CGeoPolyhedron<Real>& aNewPolyhedron, CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoPlane<Real>& aPlane, const Real aTolerance)
    {
        CGeoPlane<Real> aPlaneInv = aPlane;
        aPlaneInv.normal() = -aPlane.normal();
        aPlaneInv.setD(aPlane.d());

        aNewPolyhedron = this->clip(aNewPolygon, aNewPolygonId, aPlane, aTolerance);
        return this->clip(aNewPolygon, aNewPolygonId, aPlane, aTolerance);
    }

    template <typename Real>
    CGeoPolyhedron<Real> CGeoPolyhedron<Real>::cut(CGeoPolyhedron<Real>& aNewPolyhedron, CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoNormal<Real>& aNormal, Real& d, Real volumeFractionReq, Real& volumeFractionAct, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance)
    {
        CGeoNormal<Real> aNormalInv = -aNormal;
        aNewPolyhedron = this->clip(aNewPolygon, aNewPolygonId, aNormalInv, d, 1.0 - volumeFractionReq, volumeFractionAct, nIterations, nMaxIterations, aTolerance);
        return this->clip(aNewPolygon, aNewPolygonId, aNormal, d, volumeFractionReq, volumeFractionAct, nIterations, nMaxIterations, aTolerance);
    }

    template <typename Real>
    void CGeoPolyhedron<Real>::set(CGeoPolyhedron<Real>& aPolyhedron)
    {
        this->reset();

        for (Integer i = 0; i < aPolyhedron.nbPolygons(); ++i) {
            Integer aPolygonId = aPolyhedron.polygonId(i);

            this->addPolygon(aPolygonId, aPolyhedron.polygon(aPolygonId));
        }
    }
}
}
