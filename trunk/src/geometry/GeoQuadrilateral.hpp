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
#include "GeoTriangle.hpp"
#include "GeoVertexList.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    class CGeoQuadrilateral : public CGeoArea<Real>, public CGeoVertexList<Real> {
    public:
        CGeoQuadrilateral();
        virtual ~CGeoQuadrilateral();

        void reset();

        inline void calculateCentroid(bool bReCalculate = false);
        inline void calculateNormal(bool bReCalculate = false);
        inline void calculateArea(bool bReCalculate = false);
        inline void calculateBoundingBox(bool bReCalculate = false);

        inline bool contains(const CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance = 0.0);

        void decimate(std::vector<CGeoTriangle<Real>>& sTriangles);
    };
}
}

#include "GeoQuadrilateral_Imp.hpp"
