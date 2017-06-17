// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoLength.hpp"
#include "GeoVertexList.hpp"
#include "GeoLineList.hpp"
#include "GeoPlane.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoPolyline : public CGeoLength<Real>, public CGeoVertexList<Real> 
        {    
        public:
            CGeoPolyline();
            CGeoPolyline(CGeoLineList<Real>& aLineList, bool sort = false, const Real aTolerance = 0);

            ~CGeoPolyline();

            void set(CGeoLineList<Real>& aLineList, bool sort = false, const Real aTolerance = 0);

            Integer nbLines();
            CGeoLine<Real> line(Integer aLineIndex);

            bool isClosed(const Real aTolerance = 0);
            void close(const Real aTolerance = 0);

            inline void calculateLength(bool bReCalculate = false);
            inline void calculateBoundingBox(bool bReCalculate = false);

            inline CGeoLineList<Real> clip(CGeoPlane<Real>& aPlane);

            inline bool intersects(CGeoPlane<Real>& aPlane);

            void invert();

        };

    }

}

#include "GeoPolyline_Imp.hpp"

