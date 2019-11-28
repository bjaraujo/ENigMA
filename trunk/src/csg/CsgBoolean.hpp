// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

// Based on Evan Wallace CSG code.

#pragma once

#include "GeoPolygon.hpp"
#include "StlUtils.hpp"
#include "CsgBspNode.hpp"

namespace ENigMA
{
    namespace csg
    {
        template<typename Real> 
        class CCsgBoolean
        {
        private:

            Real m_tolerance;
            std::vector<CGeoPolygon<Real> > m_polygons;

        public:

            CCsgBoolean(const Real aTolerance = 0.0);
            explicit CCsgBoolean(std::vector<ENigMA::geometry::CGeoPolygon<Real> >& sPolygons, const Real aTolerance = 0.0);
            
            ~CCsgBoolean();

            void setTolerance(const Real aTolerance);

            Integer nbPolygons() const;
            std::vector<ENigMA::geometry::CGeoPolygon<Real> >& polygons();

            CCsgBoolean<Real> add(CCsgBoolean<Real> aCSg);
            CCsgBoolean<Real> subtract(CCsgBoolean<Real> aCSg);
            CCsgBoolean<Real> intersect(CCsgBoolean<Real> aCSg);
            
            void fromSTL(ENigMA::stl::CStlUtils<Real>& aStl);
            ENigMA::stl::CStlUtils<Real> toStl();
        };
    }
}
    
#include "CsgBoolean_Imp.hpp"
