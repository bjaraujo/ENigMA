// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "GeoBoundingBox.hpp"
#include "GeoLine.hpp"
#include "GeoVector.hpp"

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        CGeoTriangle<Real>::CGeoTriangle()
        {
        }

        template <typename Real>
        CGeoTriangle<Real>::~CGeoTriangle()
        {
        }

        template <typename Real>
        void CGeoTriangle<Real>::reset()
        {
            CGeoVertexList<Real>::reset();

            this->m_bCentroid = false;
            this->m_bNormal = false;
            this->m_bArea = false;
            this->m_bBoundingBox = false;
        }

        template <typename Real>
        CGeoPlane<Real> CGeoTriangle<Real>::getPlane()
        {
            Real d = this->normal().dot(this->vertex(0));

            CGeoPlane<Real> aPlane(this->normal(), d);

            return aPlane;
        }

        template <typename Real>
        void CGeoTriangle<Real>::calculateCentroid(bool bReCalculate)
        {
            if (!this->m_bCentroid || bReCalculate)
            {
                this->m_centroid = this->m_vertices[0];
                this->m_centroid += this->m_vertices[1];
                this->m_centroid += this->m_vertices[2];

                this->m_centroid /= 3.0;

                this->m_bCentroid = true;
            }
        }

        template <typename Real>
        void CGeoTriangle<Real>::calculateNormal(bool bReCalculate)
        {
            if (!this->m_bNormal || bReCalculate)
            {
                this->m_normal = (this->m_vertices[1] - this->m_vertices[0]).cross(this->m_vertices[2] - this->m_vertices[0]);

                this->m_area = this->m_normal.norm() * 0.5;

                if (this->m_area > 0.0)
                    this->m_normal.normalize();

                this->m_bNormal = true;
                this->m_bArea = true;
            }
        }

        template <typename Real>
        void CGeoTriangle<Real>::calculateArea(bool bReCalculate)
        {
            this->calculateNormal(bReCalculate);
        }

        template <typename Real>
        void CGeoTriangle<Real>::calculateBoundingBox(bool bReCalculate)
        {
            if (!this->m_bBoundingBox || bReCalculate)
            {
                this->m_boundingBox.reset();

                for (Integer k = 0; k < 3; ++k)
                    this->m_boundingBox.addCoordinate(this->m_vertices[k]);

                this->m_bBoundingBox = true;
            }
        }

        template <typename Real>
        bool CGeoTriangle<Real>::intersects(CGeoTriangle<Real>& aTriangle, CGeoIntersectionType& anIntersectionType, const Real aTolerance)
        {
            anIntersectionType = IT_NONE;

            // Check if bounding boxes intersect
            this->calculateBoundingBox();
            aTriangle.calculateBoundingBox();

            if (!this->boundingBox().intersects(aTriangle.boundingBox(), aTolerance))
                return false;

            this->calculateNormal();

            Real d1 = this->normal().dot(this->vertex(0));

            Real s1[3];

            for (char i = 0; i < 3; ++i)
                s1[i] = this->normal().dot(aTriangle.vertex(i)) - d1;

            if ((s1[0] > +aTolerance && s1[1] > +aTolerance && s1[2] > +aTolerance) || (s1[0] < -aTolerance && s1[1] < -aTolerance && s1[2] < -aTolerance))
            {
                return false;
            }

            aTriangle.calculateNormal();

            Real d2 = aTriangle.normal().dot(aTriangle.vertex(0));

            Real s2[3];

            for (char i = 0; i < 3; ++i)
                s2[i] = aTriangle.normal().dot(aTriangle.vertex(i)) - d2;

            if ((s2[0] > +aTolerance && s2[1] > +aTolerance && s2[2] > +aTolerance) || (s2[0] < -aTolerance && s2[1] < -aTolerance && s2[2] < -aTolerance))
            {
                return false;
            }

            Integer nVertex = 0;

            for (Integer i = 0; i < 3; ++i)
            {
                for (Integer j = 0; j < 3; ++j)
                {
                    Real d = (this->vertex(i) - aTriangle.vertex(j)).norm();

                    if (d < aTolerance)
                        nVertex++;
                }
            }

            if (nVertex == 3)
            {
                anIntersectionType = IT_COINCIDENT;
                return true;
            }

            Real c = this->normal().cross(aTriangle.normal()).norm();

            if (c <= aTolerance && fabs(d1 - d2) <= aTolerance)
            {
                // Triangles are co-planar
                for (Integer k = 0; k < 3; ++k)
                {
                    CGeoLine<Real> aLine1(this->vertex((k + 0) % 3), this->vertex((k + 1) % 3));

                    for (Integer l = 0; l < 3; l++)
                    {
                        CGeoLine<Real> aLine2(aTriangle.vertex((l + 0) % 3), aTriangle.vertex((l + 1) % 3));

                        CGeoIntersectionType aLineIntersectionType;

                        if (aLine1.intersects(aLine2, aLineIntersectionType, aTolerance))
                        {
                            if (aLineIntersectionType == IT_INTERNAL)
                            {
                                if (nVertex == 2)
                                    anIntersectionType = IT_SWAP;
                                else
                                    anIntersectionType = IT_INTERNAL;

                                return true;
                            }
                        }
                    }
                }

                if (nVertex == 2)
                {
                    anIntersectionType = IT_EDGE;
                    return true;
                }
            }
            else
            {
                if (nVertex == 2)
                {
                    anIntersectionType = IT_EDGE;
                    return true;
                }

                CGeoCoordinate<Real> aPoint;

                CGeoPlane<Real> aPlane1(this->normal(), d1);

                for (Integer k = 0; k < 3; ++k)
                {
                    CGeoLine<Real> aLine(aTriangle.vertex((k + 0) % 3), aTriangle.vertex((k + 1) % 3));

                    CGeoIntersectionType aLineIntersectionType1;

                    if (aLine.intersects(aPlane1, aPoint, aLineIntersectionType1, aTolerance))
                    {
                        CGeoIntersectionType aLineIntersectionType2(IT_NONE);

                        if (this->contains(aPoint, aLineIntersectionType2, aTolerance))
                        {
                            if (aLineIntersectionType2 == IT_INTERNAL || aLineIntersectionType2 == IT_EDGE)
                            {
                                anIntersectionType = IT_INTERNAL;
                                return true;
                            }
                        }
                    }
                }

                CGeoPlane<Real> aPlane2(aTriangle.normal(), d2);

                for (Integer k = 0; k < 3; ++k)
                {
                    CGeoLine<Real> aLine(this->vertex((k + 0) % 3), this->vertex((k + 1) % 3));

                    CGeoIntersectionType aLineIntersectionType1;

                    if (aLine.intersects(aPlane2, aPoint, aLineIntersectionType1, aTolerance))
                    {
                        CGeoIntersectionType aLineIntersectionType2;

                        if (aTriangle.contains(aPoint, aLineIntersectionType2, aTolerance))
                        {
                            if (aLineIntersectionType2 == IT_INTERNAL || aLineIntersectionType2 == IT_EDGE)
                            {
                                anIntersectionType = IT_INTERNAL;
                                return true;
                            }
                        }
                    }
                }
            }

            if (nVertex == 1)
            {
                anIntersectionType = IT_VERTEX;
                return true;
            }

            CGeoIntersectionType aPointIntersectionType;

            if ((this->contains(aTriangle.vertex(0), aPointIntersectionType, aTolerance) && this->contains(aTriangle.vertex(1), aPointIntersectionType, aTolerance) && this->contains(aTriangle.vertex(2), aPointIntersectionType, aTolerance)) || (aTriangle.contains(this->vertex(0), aPointIntersectionType, aTolerance) && aTriangle.contains(this->vertex(1), aPointIntersectionType, aTolerance) && aTriangle.contains(this->vertex(2), aPointIntersectionType, aTolerance)))
            {
                anIntersectionType = IT_INTERNAL;
                return true;
            }

            return false;
        }

        template <typename Real>
        bool CGeoTriangle<Real>::intersects(CGeoLine<Real>& aLine, CGeoCoordinate<Real>& aNewPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance)
        {
            CGeoPlane<Real> aPlane(this->normal(), this->normal().dot(this->vertex(0)));

            if (aLine.intersects(aPlane, aNewPoint, anIntersectionType, aTolerance))
            {
                if (this->contains(aNewPoint, anIntersectionType, aTolerance))
                    return true;
            }

            return false;
        }

        template <typename Real>
        bool CGeoTriangle<Real>::contains(const CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance)
        {
            anIntersectionType = IT_NONE;

            this->calculateBoundingBox();

            if (!this->boundingBox().contains(aPoint, aTolerance))
                return false;

            this->calculateNormal();

            Integer w = 0;

            Real smin = std::numeric_limits<Real>::max();

            Real s[3];

            for (Integer k = 0; k < 3; ++k)
            {
                CGeoVector<Real> v = aPoint - this->m_vertices[(k + 0) % 3];

                if (v.norm() <= aTolerance)
                {
                    anIntersectionType = IT_VERTEX;
                    return true;
                }

                s[k] = (this->m_vertices[(k + 1) % 3] - this->m_vertices[(k + 0) % 3]).cross(v).dot(this->normal());

                if (s[k] <= -aTolerance)
                    return false;

                if (s[k] < smin)
                {
                    smin = s[k];
                    w = k;
                }
            }

            if (fabs(s[(w + 0) % 3] - (s[(w + 1) % 3] + s[(w + 2) % 3])) < aTolerance)
                anIntersectionType = IT_EDGE;
            else
                anIntersectionType = IT_INTERNAL;

            return true;
        }

        template <typename Real>
        bool CGeoTriangle<Real>::distance(const CGeoCoordinate<Real>& aPoint, CGeoCoordinate<Real>& aNewPoint, Real& aDistance, const Real aTolerance)
        {
            this->calculateNormal();

            CGeoPlane<Real> aPlane(this->normal(), this->normal().dot(this->vertex(0)));

            aPlane.distance(aPoint, aNewPoint, aDistance, aTolerance);

            CGeoIntersectionType anIntersectionType;

            if (!this->contains(aNewPoint, anIntersectionType, aTolerance))
            {
                aDistance = std::numeric_limits<Real>::max();

                // Get shortest distance
                for (Integer i = 0; i < 3; ++i)
                {
                    CGeoCoordinate<Real> aVertex = this->vertex(i);

                    Real dist = (aPoint - aVertex).norm();

                    if (dist < aDistance)
                    {
                        aDistance = dist;
                        aNewPoint = aVertex;
                    }
                }
            }

            return true;
        }
    }
}
