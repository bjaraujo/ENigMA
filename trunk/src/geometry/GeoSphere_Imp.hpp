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

    namespace geometry
    {

        template <typename Real>
        CGeoSphere<Real>::CGeoSphere(const CGeoCoordinate<Real>& aPoint1, const CGeoCoordinate<Real>& aPoint2, const CGeoCoordinate<Real>& aPoint3, const CGeoCoordinate<Real>& aPoint4, const Real aTolerance)
        {

            std::vector<CGeoCoordinate<Real> > sPoints;

            sPoints.push_back(aPoint1);
            sPoints.push_back(aPoint2);
            sPoints.push_back(aPoint3);
            sPoints.push_back(aPoint4);

            Eigen::Matrix<Real, 4, 4> a;

            for (Integer i = 0; i < 4; ++i)
                a.row(i) << sPoints[i].x(), sPoints[i].y(), sPoints[i].z(), 1.0;

            Real m11 = a.determinant();

            for (Integer i = 0; i < 4; ++i)
                a.row(i) << sPoints[i].squaredNorm(), sPoints[i].y(), sPoints[i].z(), 1.0;

            Real m12 = a.determinant();

            for (Integer i = 0; i < 4; ++i)
                a.row(i) << sPoints[i].x(), sPoints[i].squaredNorm(), sPoints[i].z(), 1.0;

            Real m13 = a.determinant();

            for (Integer i = 0; i < 4; ++i)
                a.row(i) << sPoints[i].x(), sPoints[i].y(), sPoints[i].squaredNorm(), 1.0;

            Real m14 = a.determinant();

            for (Integer i = 0; i < 4; ++i)
                a.row(i) << sPoints[i].squaredNorm(), sPoints[i].x(), sPoints[i].y(), sPoints[i].z();

            Real m15 = a.determinant();

            if (fabs(m11) <= aTolerance * aTolerance)
            {
                this->m_center << 0.0, 0.0, 0.0;
                this->m_radius = 0.0;
            }
            else
            {
                this->m_center << 0.5 * m12 / m11, 0.5 * m13 / m11, 0.5 * m14 / m11;
                this->m_radius = sqrt(m_center.squaredNorm() - m15 / m11);
            }

        }

        template <typename Real>
        CGeoSphere<Real>::CGeoSphere(const CGeoCoordinate<Real>& aCenter, const Real aRadius)
        {

            this->m_center = aCenter;
            this->m_radius = aRadius;

        }
    
        template <typename Real>
        CGeoSphere<Real>::~CGeoSphere()
        {

        }

        template <typename Real>
        CGeoCoordinate<Real>& CGeoSphere<Real>::center()
        {

            return this->m_center;

        }

        template <typename Real>
        Real CGeoSphere<Real>::radius()
        {

            return this->m_radius;

        }

        template <typename Real>
        void CGeoSphere<Real>::calculateCentroid(bool bReCalculate)
        {

            if (!this->m_bCentroid || bReCalculate)
            {

                CGeoVolume<Real>::centroid() = this->m_center;

                this->m_bCentroid = true;
            }

        }

        template <typename Real>
        void CGeoSphere<Real>::calculateSurfaceArea(bool bReCalculate)
        {

            if (!this->m_bSurfaceArea || bReCalculate)
            {

                const Real pi = std::acos(-1.0);
                CGeoVolume<Real>::surfaceArea() += 4.0 * pi * this->m_radius * this->m_radius;

                this->m_bSurfaceArea = true;

            }

        }

        template <typename Real>
        void CGeoSphere<Real>::calculateVolume(bool bReCalculate)
        {

            if (!this->m_bVolume || bReCalculate)
            {

                const Real pi = std::acos(-1.0);
                CGeoVolume<Real>::volume() = 4.0 / 3.0 * pi * this->m_radius * this->m_radius * this->m_radius;

                this->m_bVolume = true;

            }

        }

        template <typename Real>
        void CGeoSphere<Real>::calculateBoundingBox(bool bReCalculate)
        {

            if (!this->m_bBoundingBox || bReCalculate)
            {

                CGeoVolume<Real>::boundingBox().reset();

                CGeoCoordinate<Real> aCoordinate1(-this->m_radius, -this->m_radius, -this->m_radius);
                CGeoCoordinate<Real> aCoordinate2(+this->m_radius, +this->m_radius, +this->m_radius);

                CGeoVolume<Real>::boundingBox().addCoordinate(aCoordinate1);
                CGeoVolume<Real>::boundingBox().addCoordinate(aCoordinate2);

                this->m_bBoundingBox = true;

            }

        }
    
        template <typename Real>
        bool CGeoSphere<Real>::contains(const CGeoCoordinate<Real>& aPoint, const Real aTolerance)
        {

            Real aDistance = (aPoint - this->m_center).norm();

            if (aDistance <= this->m_radius + aTolerance)
                return true;
            else
                return false;
        
        }

    }

}
