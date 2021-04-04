// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

namespace ENigMA {
namespace geometry {
    template <typename Real>
    CGeoCircle<Real>::CGeoCircle(const CGeoCoordinate<Real>& aPoint1, const CGeoCoordinate<Real>& aPoint2, const CGeoCoordinate<Real>& aPoint3, const Real aTolerance)
    {
        Real xDelta_a = aPoint2.x() - aPoint1.x();
        Real yDelta_a = aPoint2.y() - aPoint1.y();

        Real xDelta_b = aPoint3.x() - aPoint2.x();
        Real yDelta_b = aPoint3.y() - aPoint2.y();

        if (fabs(xDelta_a) <= aTolerance * aTolerance && fabs(yDelta_b) <= aTolerance * aTolerance) {
            this->m_center.x() = 0.5 * (aPoint2.x() + aPoint3.x());
            this->m_center.y() = 0.5 * (aPoint1.y() + aPoint2.y());
            this->m_center.z() = aPoint1.z();
            this->m_radius = (aPoint1 - this->m_center).norm();

        } else {
            Real aSlope = yDelta_a / xDelta_a; //
            Real bSlope = yDelta_b / xDelta_b;

            if (fabs(aSlope - bSlope) <= aTolerance * aTolerance) {
                this->m_center << 0.0, 0.0, 0.0;
                this->m_radius = 0.0;
            } else {
                this->m_center.x() = (aSlope * bSlope * (aPoint1.y() - aPoint3.y()) + bSlope * (aPoint1.x() + aPoint2.x()) - aSlope * (aPoint2.x() + aPoint3.x())) / (2.0 * (bSlope - aSlope));
                this->m_center.y() = -(m_center.x() - (aPoint1.x() + aPoint2.x()) * 0.5) / aSlope + (aPoint1.y() + aPoint2.y()) * 0.5;
                this->m_center.z() = aPoint1.z();

                this->m_radius = (aPoint1 - this->m_center).norm();
            }
        }
    }

    template <typename Real>
    CGeoCircle<Real>::CGeoCircle(const CGeoCoordinate<Real>& aCenter, const Real aRadius)
    {
        this->m_center = aCenter;
        this->m_radius = aRadius;
    }

    template <typename Real>
    CGeoCircle<Real>::~CGeoCircle()
    {
    }

    template <typename Real>
    const CGeoCoordinate<Real>& CGeoCircle<Real>::center() const
    {
        return this->m_center;
    }

    template <typename Real>
    Real CGeoCircle<Real>::radius() const
    {
        return this->m_radius;
    }

    template <typename Real>
    void CGeoCircle<Real>::calculateCentroid(bool bReCalculate)
    {
        CGeoArea<Real>::centroid() = this->m_center;
    }

    template <typename Real>
    void CGeoCircle<Real>::calculateNormal(bool bReCalculate)
    {
        if (!this->m_bNormal || bReCalculate) {
            CGeoArea<Real>::normal() << 1.0, 0.0, 0.0;

            this->m_bNormal = true;
        }
    }

    template <typename Real>
    void CGeoCircle<Real>::calculateArea(bool bReCalculate)
    {
        if (!this->m_bArea || bReCalculate) {
            static const Real pi = std::acos(-1.0);
            CGeoArea<Real>::area() = pi * this->m_radius * this->m_radius;

            this->m_bArea = true;
        }
    }

    template <typename Real>
    void CGeoCircle<Real>::calculateBoundingBox(bool bReCalculate)
    {
        if (!this->m_bBoundingBox || bReCalculate) {
            CGeoArea<Real>::boundingBox().reset();

            CGeoCoordinate<Real> aCoordinate1(-this->m_radius, -this->m_radius, 0.0);
            CGeoCoordinate<Real> aCoordinate2(+this->m_radius, +this->m_radius, 0.0);

            CGeoArea<Real>::boundingBox().addCoordinate(aCoordinate1);
            CGeoArea<Real>::boundingBox().addCoordinate(aCoordinate2);

            this->m_bBoundingBox = true;
        }
    }

    template <typename Real>
    bool CGeoCircle<Real>::contains(const CGeoCoordinate<Real>& aPoint, const Real aTolerance)
    {
        Real aDistance = (aPoint - this->m_center).norm();

        if (aDistance <= this->m_radius + aTolerance)
            return true;
        else
            return false;
    }
}
}
