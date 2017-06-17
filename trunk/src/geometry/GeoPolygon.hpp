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
#include "GeoPolyline.hpp"
#include "GeoPlane.hpp"

namespace ENigMA
{

    namespace geometry
    {

        // Represents a convex polygon. The vertices used to initialize 
        // a polygon must be coplanar and form a convex loop. 

        template <typename Real>
        class CGeoPolygon : public CGeoArea<Real>
        {
        protected:

            CGeoPolyline<Real> m_polyline;

        public:
            CGeoPolygon();
            ~CGeoPolygon();

            void reset();

            CGeoPolyline<Real>& polyline();
            void setPolyline(CGeoPolyline<Real>& aPolyline);

            std::vector<CGeoTriangle<Real> > triangulate();

            inline void calculateCentroid(bool bReCalculate = false);
            inline void calculateNormal(bool bReCalculate = false);
            inline void calculateArea(bool bReCalculate = false);
            inline void calculateBoundingBox(bool bReCalculate = false);

            inline CGeoPolygon<Real> clip(CGeoPlane<Real>& aPlane, const Real aTolerance = 0);

            inline bool intersects(CGeoPlane<Real>& aPlane);

            void invert();

            bool isClosed();
            void close();

        };

    }

}

#include "GeoPolygon_Imp.hpp"
