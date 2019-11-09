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
#include "GeoVector.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    class CGeoBoundingBox : public CGeoCoordinate<Real> {
    private:
        CGeoCoordinate<Real> m_center;

        CGeoVector<Real> m_min;
        CGeoVector<Real> m_max;

    public:
        CGeoBoundingBox();
        CGeoBoundingBox(const Real aMinX, const Real aMinY, const Real aMinZ, const Real aMaxX, const Real aMaxY, const Real aMaxZ);
        virtual ~CGeoBoundingBox();

        void reset();
        void addCoordinate(CGeoCoordinate<Real>& aCoordinate);

        CGeoVector<Real>& min();
        CGeoVector<Real>& max();

        void grow(const Real anAmount);
        void grow(const Real anAmountX, const Real anAmountY, const Real anAmountZ);

        void shrink(const Real anAmount);

        inline bool intersects(CGeoBoundingBox<Real>& aBoundingBox, const Real aTolerance = 0.0);
        inline bool contains(const CGeoCoordinate<Real>& aCoordinate, const Real aTolerance = 0.0);
        inline bool contains(CGeoBoundingBox<Real>& aBoundingBox, const Real aTolerance = 0.0);
    };
}
}

#include "GeoBoundingBox_Imp.hpp"
