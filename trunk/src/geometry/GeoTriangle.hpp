// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoArea.hpp"
#include "GeoLine.hpp"
#include "GeoVertexList.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    class CGeoTriangle : public CGeoArea<Real>, public CGeoVertexList<Real> {
    public:
        CGeoTriangle();
        virtual ~CGeoTriangle();

        void reset();

        CGeoPlane<Real> getPlane();

        inline void calculateCentroid(bool bReCalculate = false);
        inline void calculateNormal(bool bReCalculate = false);
        inline void calculateArea(bool bReCalculate = false);
        inline void calculateBoundingBox(bool bReCalculate = false);

        inline bool intersects(CGeoTriangle<Real>& aTriangle, CGeoIntersectionType& anIntersectionType, const Real aTolerance = 0.0);
        inline bool intersects(CGeoLine<Real>& aLine, CGeoCoordinate<Real>& aNewPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance = 0.0);
        inline bool contains(const CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance = 0.0);
        inline bool distance(const CGeoCoordinate<Real>& aPoint, CGeoCoordinate<Real>& aNewPoint, Real& aDistance, const Real aTolerance = 0.0);
    };
}
}

#include "GeoTriangle_Imp.hpp"
