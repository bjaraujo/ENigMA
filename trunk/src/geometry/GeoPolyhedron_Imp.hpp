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
            aPolyline.addVertex(aTetrahedron.vertex(1));
            aPolyline.addVertex(aTetrahedron.vertex(2));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(0, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 2
            aPolyline.addVertex(aTetrahedron.vertex(1));
            aPolyline.addVertex(aTetrahedron.vertex(0));
            aPolyline.addVertex(aTetrahedron.vertex(3));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(1, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 3
            aPolyline.addVertex(aTetrahedron.vertex(3));
            aPolyline.addVertex(aTetrahedron.vertex(0));
            aPolyline.addVertex(aTetrahedron.vertex(2));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(2, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 4
            aPolyline.addVertex(aTetrahedron.vertex(3));
            aPolyline.addVertex(aTetrahedron.vertex(2));
            aPolyline.addVertex(aTetrahedron.vertex(1));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(3, aPolygon);
            aPolyline.reset();
            aPolygon.reset();
        }

        template <typename Real>
        CGeoPolyhedron<Real>::CGeoPolyhedron(CGeoTriangularPrism<Real>& aTriangularPrism)
        {
            CGeoPolyline<Real> aPolyline;
            CGeoPolygon<Real> aPolygon;

            // face 1
            aPolyline.addVertex(aTriangularPrism.vertex(0));
            aPolyline.addVertex(aTriangularPrism.vertex(1));
            aPolyline.addVertex(aTriangularPrism.vertex(2));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(0, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 2
            aPolyline.addVertex(aTriangularPrism.vertex(5));
            aPolyline.addVertex(aTriangularPrism.vertex(4));
            aPolyline.addVertex(aTriangularPrism.vertex(3));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(1, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 3
            aPolyline.addVertex(aTriangularPrism.vertex(1));
            aPolyline.addVertex(aTriangularPrism.vertex(0));
            aPolyline.addVertex(aTriangularPrism.vertex(3));
            aPolyline.addVertex(aTriangularPrism.vertex(4));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(2, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 4
            aPolyline.addVertex(aTriangularPrism.vertex(2));
            aPolyline.addVertex(aTriangularPrism.vertex(1));
            aPolyline.addVertex(aTriangularPrism.vertex(4));
            aPolyline.addVertex(aTriangularPrism.vertex(5));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(3, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 5
            aPolyline.addVertex(aTriangularPrism.vertex(0));
            aPolyline.addVertex(aTriangularPrism.vertex(2));
            aPolyline.addVertex(aTriangularPrism.vertex(5));
            aPolyline.addVertex(aTriangularPrism.vertex(3));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(4, aPolygon);
            aPolyline.reset();
            aPolygon.reset();
        }

        template <typename Real>
        CGeoPolyhedron<Real>::CGeoPolyhedron(CGeoHexahedron<Real>& aHexahedron)
        {
            CGeoPolyline<Real> aPolyline;
            CGeoPolygon<Real> aPolygon;

            // face 1
            aPolyline.addVertex(aHexahedron.vertex(0));
            aPolyline.addVertex(aHexahedron.vertex(1));
            aPolyline.addVertex(aHexahedron.vertex(2));
            aPolyline.addVertex(aHexahedron.vertex(3));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(0, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 2
            aPolyline.addVertex(aHexahedron.vertex(7));
            aPolyline.addVertex(aHexahedron.vertex(6));
            aPolyline.addVertex(aHexahedron.vertex(5));
            aPolyline.addVertex(aHexahedron.vertex(4));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(1, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 3
            aPolyline.addVertex(aHexahedron.vertex(5));
            aPolyline.addVertex(aHexahedron.vertex(1));
            aPolyline.addVertex(aHexahedron.vertex(0));
            aPolyline.addVertex(aHexahedron.vertex(4));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(2, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 4
            aPolyline.addVertex(aHexahedron.vertex(3));
            aPolyline.addVertex(aHexahedron.vertex(2));
            aPolyline.addVertex(aHexahedron.vertex(6));
            aPolyline.addVertex(aHexahedron.vertex(7));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(3, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 5
            aPolyline.addVertex(aHexahedron.vertex(6));
            aPolyline.addVertex(aHexahedron.vertex(2));
            aPolyline.addVertex(aHexahedron.vertex(1));
            aPolyline.addVertex(aHexahedron.vertex(5));
            aPolyline.close();

            aPolygon.setPolyline(aPolyline);
            this->addPolygon(4, aPolygon);
            aPolyline.reset();
            aPolygon.reset();

            // face 6
            aPolyline.addVertex(aHexahedron.vertex(3));
            aPolyline.addVertex(aHexahedron.vertex(7));
            aPolyline.addVertex(aHexahedron.vertex(4));
            aPolyline.addVertex(aHexahedron.vertex(0));
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
        void CGeoPolyhedron<Real>::close(CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoPlane<Real>& aPlane, const Real aTolerance)
        {
            CGeoLineList<Real> aLineList;

            for (auto it = m_polygons.begin(); it != m_polygons.end(); it++)
            {
                CGeoPolyline<Real> aPolyline = it->second.polyline();

                CGeoLine<Real> aLine;
                for (Integer j = 0; j < aPolyline.nbLines(); ++j)
                {
                    aLine = aPolyline.line(j);

                    Real d1 = aPlane.distance(aLine.startPoint());
                    Real d2 = aPlane.distance(aLine.endPoint());

                    if (d1 <= aTolerance && d2 <= aTolerance)
                    {
                        aLineList.addLine(aLine);
                    }
                }
            }

            // Add new polygon
            if (aLineList.nbLines() > 2)
            {
                CGeoPolyline<Real> aPolyline(aLineList, true, aTolerance);

                aNewPolygon.setPolyline(aPolyline);

                aNewPolygon.calculateNormal(true);

                if (aPlane.normal().dot(aNewPolygon.normal()) < 0)
                    aNewPolygon.invert();

                this->addPolygon(aNewPolygonId, aNewPolygon);
            }
        }

        template <typename Real>
        void CGeoPolyhedron<Real>::calculateCentroid(bool bReCalculate)
        {
            if (!this->m_bCentroid || bReCalculate)
            {
                this->m_centroid << 0.0, 0.0, 0.0;

                for (auto it = m_polygons.begin(); it != m_polygons.end(); it++)
                {
                    it->second.calculateCentroid();
                    this->m_centroid += it->second.centroid();
                }

                if (m_polygons.size() > 0)
                    this->m_centroid /= static_cast<Real>(m_polygonIds.size());

                this->m_bCentroid = true;
            }
        }

        template <typename Real>
        void CGeoPolyhedron<Real>::calculateSurfaceArea(bool bReCalculate)
        {
            if (!this->m_bSurfaceArea || bReCalculate)
            {
                this->m_surfaceArea = 0.0;

                for (auto it = m_polygons.begin(); it != m_polygons.end(); it++)
                {
                    it->second.calculateArea();
                    this->m_surfaceArea += it->second.area();
                }

                this->m_bSurfaceArea = true;
            }
        }

        template <typename Real>
        void CGeoPolyhedron<Real>::calculateVolume(bool bReCalculate)
        {
            if (!this->m_bVolume || bReCalculate)
            {
                // see: http://en.wikipedia.org/wiki/Polyhedron#Volume
                // see: http://g3d.sourceforge.net/

                this->m_volume = 0.0;

                if (m_polygons.size() >= 4)
                {
                    CGeoPolyline<Real> aPolyline0 = m_polygons[m_polygonIds[0]].polyline();

                    if (aPolyline0.nbVertices() > 0)
                    {
                        CGeoCoordinate<Real> v0 = aPolyline0.vertex(0);

                        CGeoPolygon<Real> aPolygon;

                        for (Integer i = 1; i < static_cast<Integer>(m_polygonIds.size()); ++i)
                        {
                            aPolygon = m_polygons[m_polygonIds[i]];
                            aPolygon.calculateArea();

                            if (aPolygon.polyline().nbVertices() > 0)
                                this->m_volume += (aPolygon.polyline().vertex(0) - v0).dot(aPolygon.normal()) * aPolygon.area();
                        }
                    }
                }

                this->m_volume /= 3.0;
                this->m_bVolume = true;
            }
        }

        template <typename Real>
        void CGeoPolyhedron<Real>::calculateBoundingBox(bool bReCalculate)
        {
            if (!this->m_bBoundingBox || bReCalculate)
            {
                CGeoVolume<Real>::boundingBox().reset();

                for (auto it = m_polygons.begin(); it != m_polygons.end(); it++)
                {
                    for (Integer j = 0; j < it->second.polyline().nbVertices(); ++j)
                    {
                        CGeoVolume<Real>::boundingBox().addCoordinate(it->second.polyline().vertex(j));
                    }
                }

                this->m_bBoundingBox = true;
            }
        }

        template <typename Real>
        CGeoPolyhedron<Real> CGeoPolyhedron<Real>::clip(CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoPlane<Real>& aPlane, const Real aTolerance)
        {
            CGeoPolyhedron<Real> aPolyhedron;

            for (auto it = m_polygons.begin(); it != m_polygons.end(); it++)
            {
                CGeoPolygon<Real> aPolygon = it->second.clip(aPlane, aTolerance);

                if (aPolygon.polyline().nbVertices() > 3)
                    aPolyhedron.addPolygon(it->first, aPolygon);
            }

            aPolyhedron.close(aNewPolygon, aNewPolygonId, aPlane, aTolerance);

            return aPolyhedron;
        }

        template <typename Real>
        CGeoPolyhedron<Real> CGeoPolyhedron<Real>::clip(CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoPlane<Real>& aPlane, Real volumeFractionReq, Real& volumeFractionAct, Integer& nIterations, const Integer nMaxIterations, const Real aNormalizedTolerance, const Real aTolerance)
        {
            CGeoPolyhedron<Real> aPolyhedron;

            if (volumeFractionReq <= 0.0)
            {
                nIterations = 0;

                volumeFractionAct = 0.0;

                return aPolyhedron;
            }

            if (volumeFractionReq >= 1.0)
            {
                nIterations = 0;

                // Copy this polyhedron
                for (auto it = m_polygons.begin(); it != m_polygons.end(); it++)
                {
                    aPolyhedron.addPolygon(it->first, it->second);
                }

                volumeFractionAct = 1.0;

                return aPolyhedron;
            }

            // Find interval a, b

            this->calculateVolume();
            Real vt = CGeoVolume<Real>::volume();

            Real a = std::numeric_limits<Real>::max();
            Real b = std::numeric_limits<Real>::lowest();

            for (auto it = m_polygons.begin(); it != m_polygons.end(); it++)
            {
                for (Integer j = 0; j < it->second.polyline().nbVertices(); ++j)
                {
                    Real d = it->second.polyline().vertex(j).dot(aPlane.normal());

                    a = std::min(a, d);
                    b = std::max(b, d);
                }
            }

            // Solve using Brent's method

            Real e = a;
            Real f = std::numeric_limits<Real>::max();

            Real fa = 0.0 - volumeFractionReq;
            Real fb = 1.0 - volumeFractionReq;

            Real fc = fa;

            Real s = 0.0;
            Real fs = 0.0;

            Real vs = 0.0;

            Real tmp;

            bool mflag = true;

            nIterations = nMaxIterations;

            for (Integer i = 0; i < nMaxIterations; ++i)
            {
                if ((fb == 0.0) || (std::fabs(fa - fb) <= aTolerance))
                {
                    nIterations = i + 1;
                    break;
                }

                if ((fa != fc) && (fb != fc))
                    // Inverse quadratic interpolation
                    s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + e * fa * fb / (fc - fa) / (fc - fb);
                else
                    // Secant Rule
                    s = b - fb * (b - a) / (fb - fa);

                tmp = (3.0 * a + b) / 4.0;

                if ((!(((s > tmp) && (s < b)) ||
                       ((s < tmp) && (s > b)))) ||
                        (mflag && (std::fabs(s - b) >= (std::fabs(b - e) * 0.5))) ||
                        (!mflag && (std::fabs(s - b) >= (std::fabs(e - f) * 0.5))) ||
                        (mflag && (std::fabs(b - e) < aTolerance)) ||
                        (!mflag && (std::fabs(e - f) < aTolerance)))
                {
                    s = (a + b) * 0.5;
                    mflag = true;
                }
                else
                {
                    mflag = false;
                }

                aPlane.setD(s);
                aPolyhedron = this->clip(aNewPolygon, aNewPolygonId, aPlane, aTolerance);
                aPolyhedron.calculateVolume(true);
                vs = aPolyhedron.volume();                
                fs = vs / vt - volumeFractionReq;

                f = e;
                e = b;
                fc = fb;

                if (fa * fs < 0)
                {
                    b = s;
                    fb = fs;
                }
                else
                {
                    a = s;
                    fa = fs;
                }

                // if |f(a)| < |f(b)| then swap (a,b) end if
                if (std::fabs(fa) < std::fabs(fb))
                {
                    tmp = a;
                    a = b;
                    b = tmp;

                    tmp = fa;
                    fa = fb;
                    fb = tmp;
                }
            }

            volumeFractionAct = vs / vt;

            return aPolyhedron;
        }

        template <typename Real>
        CGeoPolyhedron<Real> CGeoPolyhedron<Real>::cut(CGeoPolyhedron<Real>& aNewPolyhedron, CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoPlane<Real>& aPlane, const Real aTolerance)
        {
            aNewPolyhedron = this->clip(aNewPolygon, aNewPolygonId, aPlane, aTolerance);
            return this->clip(aNewPolygon, aNewPolygonId, aPlane, aTolerance);
        }

        template <typename Real>
        CGeoPolyhedron<Real> CGeoPolyhedron<Real>::cut(CGeoPolyhedron<Real>& aNewPolyhedron, CGeoPolygon<Real>& aNewPolygon, const Integer aNewPolygonId, CGeoPlane<Real>& aPlane, Real volumeFractionReq, Real& volumeFractionAct, Integer& nIterations, const Integer nMaxIterations, const Real aNormalizedTolerance, const Real aTolerance)
        {
            CGeoPlane<Real> aPlaneInv;
            aPlaneInv.normal() = -aPlane.normal();
            aPlaneInv.setD(aPlane.d());

            aNewPolyhedron = this->clip(aNewPolygon, aNewPolygonId, aPlaneInv, 1.0 - volumeFractionReq, volumeFractionAct, nIterations, nMaxIterations, aNormalizedTolerance, aTolerance);
            return this->clip(aNewPolygon, aNewPolygonId, aPlane, volumeFractionReq, volumeFractionAct, nIterations, nMaxIterations, aNormalizedTolerance, aTolerance);
        }

        template <typename Real>
        void CGeoPolyhedron<Real>::set(CGeoPolyhedron<Real>& aPolyhedron)
        {
            this->reset();

            for (Integer i = 0; i < aPolyhedron.nbPolygons(); ++i)
            {
                Integer aPolygonId = aPolyhedron.polygonId(i);

                this->addPolygon(aPolygonId, aPolyhedron.polygon(aPolygonId));
            }
        }
    }
}
