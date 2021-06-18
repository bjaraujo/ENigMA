// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoCoordinate.hpp"
#include "GeoLength.hpp"
#include "GeoPlane.hpp"
#include "GeoVector.hpp"

namespace ENigMA
{
    namespace geometry
    {
        enum EIntersectionType
        {
            IT_NONE = 0,
            IT_VERTEX,
            IT_EDGE,
            IT_COINCIDENT,
            IT_INTERNAL,
            IT_SWAP
        };

        class CGeoIntersectionType
        {
        private:
            EIntersectionType m_intersectionType;

        public:
            CGeoIntersectionType()
                : m_intersectionType(IT_NONE) {};
            explicit CGeoIntersectionType(const EIntersectionType anIntersectionType)
                : m_intersectionType(anIntersectionType) {};

            void reset()
            {
                m_intersectionType = IT_NONE;
            }

            CGeoIntersectionType& operator=(const EIntersectionType& anIntersectionType)
            {
                this->m_intersectionType = anIntersectionType;
                return *this;
            }

            bool operator==(const EIntersectionType& anIntersectionType) const
            {
                return this->m_intersectionType == anIntersectionType;
            }
        };

        template <typename Real>
        class CGeoLine : public CGeoLength<Real>
        {
        protected:
            CGeoCoordinate<Real> m_startPoint;
            CGeoCoordinate<Real> m_endPoint;

            CGeoVector<Real> m_vector;

        public:
            CGeoLine();
            CGeoLine(const CGeoCoordinate<Real>& aPoint1, const CGeoCoordinate<Real>& aPoint2);
            virtual ~CGeoLine();

            void reset();

            void setStartPoint(const CGeoCoordinate<Real>& aPoint1);
            void setEndPoint(const CGeoCoordinate<Real>& aPoint2);

            CGeoCoordinate<Real> startPoint() const;
            CGeoCoordinate<Real> endPoint() const;
            CGeoVector<Real> vector() const;

            CGeoCoordinate<Real> midPoint(Real factor);

            inline void calculateLength(bool bReCalculate = false);
            inline void calculateBoundingBox(bool bReCalculate = false);

            inline CGeoLine<Real> clip(const CGeoPlane<Real>& aPlane, const Real aTolerance = 0.0);

            inline bool intersects(const CGeoPlane<Real>& aPlane, CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance = 0.0);

            inline bool intersects(CGeoLine<Real>& aLine, CGeoIntersectionType& anIntersectionType, const Real aTolerance = 0.0);
            inline bool intersects(CGeoLine<Real>& aLine, CGeoCoordinate<Real>& aPoint, const Real aTolerance = 0.0);
            inline bool intersects(CGeoLine<Real>& aLine, CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance = 0.0);

            inline bool contains(const CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance = 0.0);

            inline bool distance(const CGeoCoordinate<Real>& aPoint, CGeoCoordinate<Real>& aNewPoint, Real& aDistance, const Real aTolerance = 0.0);
            inline bool distance(const CGeoLine<Real>& aLine, CGeoCoordinate<Real>& aPoint1, CGeoCoordinate<Real>& aPoint2, Real& aDistance, const Real aTolerance = 0.0);

            void invert();
        };
    }
}

#include "GeoLine_Imp.hpp"
