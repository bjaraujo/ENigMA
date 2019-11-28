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

namespace ENigMA
{
    namespace csg
    {
        template<typename Real> 
        class CCsgBspNode
        {
        private:

            Real m_tolerance;

            CCsgBspNode<Real>* m_front;
            CCsgBspNode<Real>* m_back;

            ENigMA::geometry::CGeoPlane<Real> m_plane;

            std::vector<ENigMA::geometry::CGeoPolygon<Real>> m_polygons;
            
        public:

            CCsgBspNode(std::vector<ENigMA::geometry::CGeoPolygon<Real>>& sPolygons, const Real aTolerance);
            explicit CCsgBspNode(const Real aTolerance);

            ~CCsgBspNode();

            void invert();

            std::vector<ENigMA::geometry::CGeoPolygon<Real>> clipPolygons(std::vector<ENigMA::geometry::CGeoPolygon<Real>>& sPolygons, bool removeCoplanarFront = false);

            void clipTo(CCsgBspNode<Real> aCsgBspNode, bool removeCoplanarFront = false);

            std::vector<ENigMA::geometry::CGeoPolygon<Real>>& allPolygons();

            void build(std::vector<ENigMA::geometry::CGeoPolygon<Real>>& sPolygons);

            void splitPolygon(ENigMA::geometry::CGeoPolygon<Real> aPolygon, std::vector<ENigMA::geometry::CGeoPolygon<Real>>& sCoplanarFront, std::vector<CGeoPolygon<Real>>& sCoplanarBack, std::vector<CGeoPolygon<Real>>& sFront, std::vector<CGeoPolygon<Real>>& sBack);

            void clear();

        };
    }
}
    
#include "CsgBspNode_Imp.hpp"
