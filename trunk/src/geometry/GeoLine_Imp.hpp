// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iostream>

#include "GeoBoundingBox.hpp"

namespace ENigMA {
namespace geometry {
    template <typename Real>
    CGeoLine<Real>::CGeoLine()
    {
    }

    template <typename Real>
    CGeoLine<Real>::CGeoLine(CGeoCoordinate<Real>& aPoint1, CGeoCoordinate<Real>& aPoint2)
    {
        this->m_startPoint = aPoint1;
        this->m_endPoint = aPoint2;
        this->m_vector = aPoint2 - aPoint1;
    }

    template <typename Real>
    CGeoLine<Real>::~CGeoLine()
    {
    }

    template <typename Real>
    void CGeoLine<Real>::reset()
    {
        this->m_bLength = false;
    }

    template <typename Real>
    void CGeoLine<Real>::setStartPoint(const CGeoCoordinate<Real>& aPoint1)
    {
        this->m_startPoint = aPoint1;
    }

    template <typename Real>
    void CGeoLine<Real>::setEndPoint(const CGeoCoordinate<Real>& aPoint2)
    {
        this->m_endPoint = aPoint2;
        this->m_vector = aPoint2 - m_startPoint;
    }

    template <typename Real>
    CGeoCoordinate<Real>& CGeoLine<Real>::startPoint()
    {
        return this->m_startPoint;
    }

    template <typename Real>
    CGeoCoordinate<Real>& CGeoLine<Real>::endPoint()
    {
        return this->m_endPoint;
    }

    template <typename Real>
    CGeoVector<Real>& CGeoLine<Real>::vector()
    {
        return this->m_vector;
    }

    template <typename Real>
    CGeoCoordinate<Real> CGeoLine<Real>::midPoint(Real factor)
    {
        return this->m_startPoint + this->m_vector * factor;
    }

    template <typename Real>
    void CGeoLine<Real>::calculateLength(bool bReCalculate)
    {
        if (!this->m_bLength || bReCalculate) {
            this->m_length = m_vector.norm();

            this->m_bLength = true;
        }
    }

    template <typename Real>
    void CGeoLine<Real>::calculateBoundingBox(bool bReCalculate)
    {
        if (!this->m_bBoundingBox || bReCalculate) {
            this->m_boundingBox.reset();

            this->m_boundingBox.addCoordinate(this->m_startPoint);
            this->m_boundingBox.addCoordinate(this->m_endPoint);

            this->m_bBoundingBox = true;
        }
    }

    template <typename Real>
    CGeoLine<Real> CGeoLine<Real>::clip(CGeoPlane<Real>& aPlane)
    {
        CGeoLine<Real> aLine;
        CGeoCoordinate<Real> p1, p2;

        // Zero length line
        p1 = m_startPoint;
        p2 = m_startPoint;

        Real den = aPlane.normal().dot(m_vector);
        Real num = aPlane.normal().dot(m_startPoint);

        if (std::fabs(den) > std::numeric_limits<Real>::epsilon()) {
            Real alpha = (aPlane.d() - num) / den;

            if (alpha >= 0 && alpha <= 1) {
                if (den > 0) {
                    // Partial line
                    p1 = this->m_startPoint;
                    p2 = this->m_startPoint + this->m_vector * alpha;
                } else {
                    // Partial line
                    p1 = this->m_startPoint + this->m_vector * alpha;
                    p2 = this->m_startPoint + this->m_vector;
                }

            } else {
                if (num < aPlane.d()) {
                    // Original line
                    p1 = this->m_startPoint;
                    p2 = this->m_startPoint + this->m_vector;
                }
            }

        } else {
            if (num < aPlane.d()) {
                // Original line
                p1 = this->m_startPoint;
                p2 = this->m_startPoint + this->m_vector;
            }
        }

        aLine.setStartPoint(p1);
        aLine.setEndPoint(p2);

        return aLine;
    }

    template <typename Real>
    bool CGeoLine<Real>::intersects(CGeoPlane<Real>& aPlane, CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance)
    {
        anIntersectionType = IT_NONE;

        CGeoCoordinate<Real> p3 = aPlane.normal() * aPlane.d();

        CGeoVector<Real> u = p3 - this->m_startPoint;
        CGeoVector<Real> v = this->m_vector;

        Real d = aPlane.normal().dot(v);

        if (fabs(d) <= aTolerance * aTolerance)
            return false;

        Real s = aPlane.normal().dot(u) / d;

        if (s >= -aTolerance && s <= 1.0 + aTolerance) {
            aPoint = this->m_startPoint + s * this->m_vector;

            if (s <= aTolerance || s >= 1.0 - aTolerance)
                anIntersectionType = IT_VERTEX;
            else
                anIntersectionType = IT_INTERNAL;

            return true;
        }

        return false;
    }

    template <typename Real>
    bool CGeoLine<Real>::intersects(CGeoLine<Real>& aLine, CGeoIntersectionType& anIntersectionType, const Real aTolerance)
    {
        CGeoCoordinate<Real> aPoint;

        return this->intersects(aLine, aPoint, anIntersectionType, aTolerance);
    }

    template <typename Real>
    bool CGeoLine<Real>::intersects(CGeoLine<Real>& aLine, CGeoCoordinate<Real>& aPoint, const Real aTolerance)
    {
        CGeoIntersectionType anIntersectionType;

        return this->intersects(aLine, aPoint, anIntersectionType, aTolerance);
    }

    template <typename Real>
    bool CGeoLine<Real>::intersects(CGeoLine<Real>& aLine, CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance)
    {
        // http://mathworld.wolfram.com/Line-LineIntersection.html
        // in 3d; will also work in 2d if z components are 0

        anIntersectionType = IT_NONE;

        // Check if bounding boxes intersect
        this->calculateBoundingBox();
        aLine.calculateBoundingBox();

        if (!this->boundingBox().intersects(aLine.boundingBox(), aTolerance))
            return false;

        CGeoVector<Real>& da = this->m_vector;
        CGeoVector<Real>& db = aLine.vector();

        CGeoVector<Real> dc = aLine.startPoint() - this->m_startPoint;
        CGeoVector<Real> dd = da.cross(db);

        Real a = fabs(dc.dot(dd));

        if (a > aTolerance) // lines are not coplanar
            return false;

        Real d1 = (this->startPoint() - aLine.startPoint()).norm();
        Real d2 = (this->endPoint() - aLine.endPoint()).norm();

        if (d1 <= aTolerance && d2 <= aTolerance) {
            aPoint = this->m_startPoint;
            anIntersectionType = IT_COINCIDENT;
            return true;
        }

        Real d3 = (this->startPoint() - aLine.endPoint()).norm();
        Real d4 = (this->endPoint() - aLine.startPoint()).norm();

        if (d3 <= aTolerance && d4 <= aTolerance) {
            aPoint = this->m_startPoint;
            anIntersectionType = IT_COINCIDENT;
            return true;
        }

        Real d = dd.squaredNorm();

        if (d <= aTolerance * aTolerance) {
            if (d1 <= aTolerance || d2 <= aTolerance) {
                aPoint = this->m_startPoint;

                if (+db.dot(da) >= aTolerance * aTolerance) {
                    anIntersectionType = IT_INTERNAL;
                    return true;
                } else {
                    if (d2 <= aTolerance)
                        aPoint = this->m_endPoint;

                    anIntersectionType = IT_VERTEX;
                    return true;
                }
            }

            if (d3 <= aTolerance || d4 <= aTolerance) {
                aPoint = this->m_startPoint;

                if (-db.dot(da) >= aTolerance * aTolerance) {
                    anIntersectionType = IT_INTERNAL;
                    return true;
                } else {
                    if (d4 <= aTolerance)
                        aPoint = this->m_endPoint;

                    anIntersectionType = IT_VERTEX;
                    return true;
                }
            }

            if (this->contains(aLine.startPoint(), anIntersectionType, aTolerance)) {
                aPoint = aLine.startPoint();
                anIntersectionType = IT_INTERNAL;
                return true;
            }

            if (this->contains(aLine.endPoint(), anIntersectionType, aTolerance)) {
                aPoint = aLine.endPoint();
                anIntersectionType = IT_INTERNAL;
                return true;
            }

            if (aLine.contains(this->startPoint(), anIntersectionType, aTolerance)) {
                aPoint = this->startPoint();
                anIntersectionType = IT_INTERNAL;
                return true;
            }

            if (aLine.contains(this->endPoint(), anIntersectionType, aTolerance)) {
                aPoint = this->endPoint();
                anIntersectionType = IT_INTERNAL;
                return true;
            }

        } else {
            if (d1 <= aTolerance || d2 <= aTolerance) {
                aPoint = this->m_startPoint;

                if (d2 <= aTolerance)
                    aPoint = this->m_endPoint;

                anIntersectionType = IT_VERTEX;
                return true;
            }

            if (d3 <= aTolerance || d4 <= aTolerance) {
                aPoint = this->m_startPoint;

                if (d4 <= aTolerance)
                    aPoint = this->m_endPoint;

                anIntersectionType = IT_VERTEX;
                return true;
            }

            Real s = dc.cross(db).dot(dd) / d;

            if (s >= -aTolerance && s <= 1.0 + aTolerance) {
                aPoint = this->m_startPoint + s * this->m_vector;

                if (aLine.contains(aPoint, anIntersectionType, aTolerance)) {
                    if (s > aTolerance && s < 1.0 - aTolerance)
                        anIntersectionType = IT_INTERNAL;

                    return true;
                }
            }
        }

        return false;
    }

    template <typename Real>
    bool CGeoLine<Real>::contains(const CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance)
    {
        anIntersectionType = IT_NONE;

        CGeoVector<Real> v = aPoint - m_startPoint;

        Real angle = m_vector.angle(v);

        if (angle > 3.14 / 180.0) // 1 degree
            return false;

        CGeoNormal<Real> n = m_vector;
        n.normalize();

        Real s = n.dot(v) / m_vector.norm();

        if (s >= -aTolerance && s <= 1.0 + aTolerance) {
            if (s <= aTolerance || s >= 1.0 - aTolerance)
                anIntersectionType = IT_VERTEX;
            else
                anIntersectionType = IT_INTERNAL;

            return true;
        }

        return false;
    }

    template <typename Real>
    inline bool CGeoLine<Real>::distance(const CGeoCoordinate<Real>& aPoint, CGeoCoordinate<Real>& aNewPoint, Real& aDistance, const Real aTolerance)
    {
        CGeoVector<Real> v = this->m_vector;
        CGeoVector<Real> w = aPoint - this->startPoint();

        Real d = v.dot(v);

        if (fabs(d) > aTolerance * aTolerance) {
            Real n = w.dot(v);
            Real s = n / d;

            if (s >= -aTolerance && s <= 1.0 + aTolerance)
                aNewPoint = this->startPoint() + s * v;
            else if (s < -aTolerance)
                aNewPoint = this->startPoint();
            else if (s > 1.0 + aTolerance)
                aNewPoint = this->endPoint();

            aDistance = (aPoint - aNewPoint).norm();

            return true;
        }

        return false;
    }

    template <typename Real>
    inline bool CGeoLine<Real>::distance(CGeoLine<Real>& aLine, CGeoCoordinate<Real>& aPoint1, CGeoCoordinate<Real>& aPoint2, Real& aDistance, const Real aTolerance)
    {
        // http://paulbourke.net/geometry/pointlineplane/

        CGeoVector<Real>& v21 = this->vector();
        CGeoVector<Real>& v43 = aLine.vector();

        Real d4343 = v43.dot(v43);

        CGeoVector<Real> v13 = this->startPoint() - aLine.startPoint();

        Real d1343 = v13.dot(v43);
        Real d4321 = v43.dot(v21);
        Real d1321 = v13.dot(v21);
        Real d2121 = v21.dot(v21);

        Real den = d2121 * d4343 - d4321 * d4321;
        Real num = d1343 * d4321 - d1321 * d4343;

        if (fabs(den) <= aTolerance)
            return false;

        Real mua = num / den;

        if (mua < -aTolerance)
            mua = 0.0;

        if (mua > 1.0 + aTolerance)
            mua = 1.0;

        Real mub = (d1343 + d4321 * mua) / d4343;

        aPoint1 = this->startPoint() + mua * this->vector();
        aPoint2 = aLine.startPoint() + mub * aLine.vector();

        aDistance = (aPoint2 - aPoint1).norm();

        return true;
    }

    template <typename Real>
    void CGeoLine<Real>::invert()
    {
        m_startPoint += m_vector;
        m_endPoint -= m_vector;

        m_vector *= -1;
    }
}
}
