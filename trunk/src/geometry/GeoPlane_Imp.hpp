// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

// <ProjectName> $ProjectName$ </ProjectName>

#pragma once

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        CGeoPlane<Real>::CGeoPlane()
        {
        }

        template <typename Real>
        CGeoPlane<Real>::CGeoPlane(const CGeoNormal<Real>& aNormal, const Real d)
        {
            this->set(aNormal, d);
        }

        template <typename Real>
        void CGeoPlane<Real>::set(const CGeoNormal<Real>& aNormal, const Real d)
        {
            m_normal = aNormal;
            m_d = d;
        }

        template <typename Real>
        CGeoPlane<Real>::~CGeoPlane()
        {
        }

        template <typename Real>
        CGeoNormal<Real>& CGeoPlane<Real>::normal()
        {
            return m_normal;
        }

        template <typename Real>
        void CGeoPlane<Real>::setD(const Real aValue)
        {
            m_d = aValue;
        }

        template <typename Real>
        Real CGeoPlane<Real>::d() const
        {
            return m_d;
        }

        template <typename Real>
        void CGeoPlane<Real>::invert()
        {
            m_normal *= -1;
            m_d *= -1;
        };

        template <typename Real>
        Real CGeoPlane<Real>::distance(const CGeoCoordinate<Real>& aPoint)
        {
            Real s = this->normal().dot(aPoint) - this->d();
            return fabs(s);
        }

        template <typename Real>
        Real CGeoPlane<Real>::distance(const CGeoCoordinate<Real>& aPoint, CGeoCoordinate<Real>& aNewPoint)
        {
            Real s = this->normal().dot(aPoint) - this->d();
            aNewPoint = aPoint - this->normal() * s;
            return fabs(s);
        }
    }
}
