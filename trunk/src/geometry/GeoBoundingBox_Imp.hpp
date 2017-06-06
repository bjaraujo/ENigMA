// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#pragma once

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        CGeoBoundingBox<Real>::CGeoBoundingBox()
        {

            this->reset();

        }

        template <typename Real>
        CGeoBoundingBox<Real>::CGeoBoundingBox(const Real aMinX, const Real aMinY, const Real aMinZ, const Real aMaxX, const Real aMaxY, const Real aMaxZ)
        {

            m_min.x() = aMinX;
            m_min.y() = aMinY;
            m_min.z() = aMinZ;

            m_max.x() = aMaxX;
            m_max.y() = aMaxY;
            m_max.z() = aMaxZ;

            // Center
            m_center = (m_min + m_max) * 0.5;

        }

        template <typename Real>
        CGeoBoundingBox<Real>::~CGeoBoundingBox()
        {

        }

        template <typename Real>
        void CGeoBoundingBox<Real>::reset()
        {

            m_min.x() = +std::numeric_limits<Real>::max();
            m_min.y() = +std::numeric_limits<Real>::max();
            m_min.z() = +std::numeric_limits<Real>::max();

            m_max.x() = -std::numeric_limits<Real>::max();
            m_max.y() = -std::numeric_limits<Real>::max();
            m_max.z() = -std::numeric_limits<Real>::max();

            m_center << 0.0, 0.0, 0.0;

        }

        template <typename Real>
        void CGeoBoundingBox<Real>::addCoordinate(CGeoCoordinate<Real>& aCoordinate)
        {

            m_min.x() = std::min(m_min.x(), aCoordinate.x());
            m_min.y() = std::min(m_min.y(), aCoordinate.y());
            m_min.z() = std::min(m_min.z(), aCoordinate.z());

            m_max.x() = std::max(m_max.x(), aCoordinate.x());
            m_max.y() = std::max(m_max.y(), aCoordinate.y());
            m_max.z() = std::max(m_max.z(), aCoordinate.z());

            // Center
            m_center = (m_min + m_max) * 0.5;

        }

        template <typename Real>
        CGeoVector<Real>& CGeoBoundingBox<Real>::min()
        {

            return m_min;

        }

        template <typename Real>
        CGeoVector<Real>& CGeoBoundingBox<Real>::max()
        {

            return m_max;

        }

        template <typename Real>
        bool CGeoBoundingBox<Real>::intersects(CGeoBoundingBox<Real>& aBoundingBox, const Real aTolerance)
        {

            //http://www.miguelcasillas.com/?p=30

            return (aBoundingBox.max().x() >= this->m_min.x() - aTolerance && 
                    aBoundingBox.min().x() <= this->m_max.x() + aTolerance &&
                    aBoundingBox.max().y() >= this->m_min.y() - aTolerance &&
                    aBoundingBox.min().y() <= this->m_max.y() + aTolerance &&
                    aBoundingBox.max().z() >= this->m_min.z() - aTolerance &&
                    aBoundingBox.min().z() <= this->m_max.z() + aTolerance);

        }

        template <typename Real>
        bool CGeoBoundingBox<Real>::contains(const CGeoCoordinate<Real>& aCoordinate, const Real aTolerance)
        {

            return (aCoordinate.x() >= this->m_min.x() - aTolerance && 
                    aCoordinate.x() <= this->m_max.x() + aTolerance &&
                    aCoordinate.y() >= this->m_min.y() - aTolerance &&
                    aCoordinate.y() <= this->m_max.y() + aTolerance &&
                    aCoordinate.z() >= this->m_min.z() - aTolerance &&
                    aCoordinate.z() <= this->m_max.z() + aTolerance);

        }

        template <typename Real>
        bool CGeoBoundingBox<Real>::contains(CGeoBoundingBox<Real>& aBoundingBox, const Real aTolerance)
        {

            return (this->m_min.x() - aTolerance <= aBoundingBox.min().x() && 
                    this->m_max.x() + aTolerance >= aBoundingBox.max().x() &&
                    this->m_min.y() - aTolerance <= aBoundingBox.min().y() &&
                    this->m_max.y() + aTolerance >= aBoundingBox.max().y() &&
                    this->m_min.z() - aTolerance <= aBoundingBox.min().z() &&
                    this->m_max.z() + aTolerance >= aBoundingBox.max().z());

        }

        template <typename Real>
        void CGeoBoundingBox<Real>::grow(const Real anAmount)
        {

            m_min.array() -= anAmount;
            m_max.array() += anAmount;

        }

        template <typename Real>
        void CGeoBoundingBox<Real>::shrink(const Real anAmount)
        {

            m_min.array() += anAmount;
            m_max.array() -= anAmount;

        }

    }

}
