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

using namespace ENigMA::geometry;

namespace ENigMA
{
    namespace csg
    {
        template<typename Real> 
        CCsgBoolean<Real>::CCsgBoolean(const Real aTolerance)
        {
            m_tolerance = aTolerance;
        }

        template<typename Real> 
        CCsgBoolean<Real>::CCsgBoolean(std::vector<CGeoPolygon<Real>>& sPolygons, const Real aTolerance)
        {
            m_tolerance = aTolerance;

            for (Integer i = 0; i < static_cast<Integer> (sPolygons.size()); ++i)
                m_polygons.push_back(sPolygons[i]);
        }

        template<typename Real> 
        CCsgBoolean<Real>::~CCsgBoolean()
        {
        
        }

        template<typename Real> 
        void CCsgBoolean<Real>::setTolerance(const Real aTolerance)
        {
        
            m_tolerance = aTolerance;

        }

        template<typename Real> 
        Integer CCsgBoolean<Real>::nbPolygons() const
        {
            return static_cast<Integer>(m_polygons.size());
        }
        
        template<typename Real> 
        std::vector<CGeoPolygon<Real>>& CCsgBoolean<Real>::polygons()
        {
            return this->m_polygons;
        }

        // Return a new CSG solid representing space in either this solid or in the
        // solid 'csg'. Neither this solid nor the solid 'csg' are modified.
        //
        //     A.add(B)
        //
        //     +-------+            +-------+
        //     |       |            |       |
        //     |   A   |            |       |
        //     |    +--+----+   =   |       +----+
        //     +----+--+    |       +----+       |
        //          |   B   |            |       |
        //          |       |            |       |
        //          +-------+            +-------+
        //

        template<typename Real> 
        CCsgBoolean<Real> CCsgBoolean<Real>::add(CCsgBoolean<Real> aCSg)
        {
            CCsgBspNode<Real> a(this->m_polygons, m_tolerance);
            CCsgBspNode<Real> b(aCSg.m_polygons, m_tolerance);

            a.clipTo(b);
            b.clipTo(a);
            b.invert();
            b.clipTo(a);
            b.invert();
            a.build(b.allPolygons());

            CCsgBoolean<Real> result(a.allPolygons(), m_tolerance);
            a.clear();
            b.clear();

            return result;
        };

        // Return a new CSG solid representing space in this solid but not in the
        // solid 'csg'. Neither this solid nor the solid 'csg' are modified.
        //
        //     A.subtract(B)
        //
        //     +-------+            +-------+
        //     |       |            |       |
        //     |   A   |            |       |
        //     |    +--+----+   =   |    +--+
        //     +----+--+    |       +----+
        //          |   B   |
        //          |       |
        //          +-------+

        template<typename Real> 
        CCsgBoolean<Real> CCsgBoolean<Real>::subtract(CCsgBoolean<Real> aCSg)
        {
            CCsgBspNode<Real> a(this->m_polygons, m_tolerance);
            CCsgBspNode<Real> b(aCSg.m_polygons, m_tolerance);

            a.invert();
            a.clipTo(b);
            b.clipTo(a, true);
            a.build(b.allPolygons());
            a.invert();

            CCsgBoolean<Real> result(a.allPolygons(), m_tolerance);
            a.clear();
            b.clear();

            return result;
        };

        // Return a new CSG solid representing space both this solid and in the
        // solid 'csg'. Neither this solid nor the solid 'csg' are modified.
        //
        //     A.intersect(B)
        //
        //     +-------+
        //     |       |
        //     |   A   |
        //     |    +--+----+   =   +--+
        //     +----+--+    |       +--+
        //          |   B   |
        //          |       |
        //          +-------+

        template<typename Real> 
        CCsgBoolean<Real> CCsgBoolean<Real>::intersect(CCsgBoolean<Real> aCSg)
        {
            CCsgBspNode<Real> a(this->m_polygons, m_tolerance);
            CCsgBspNode<Real> b(aCSg.m_polygons, m_tolerance);

            a.invert();
            b.clipTo(a);
            b.invert();
            a.clipTo(b);
            b.clipTo(a);
            a.build(b.allPolygons());
            a.invert();

            CCsgBoolean<Real> result(a.allPolygons(), m_tolerance);
            a.clear();
            b.clear();

            return result;
        };

        template<typename Real> 
        void CCsgBoolean<Real>::fromSTL(ENigMA::stl::CStlUtils<Real>& aStl)
        {
            for (Integer i = 0; i < aStl.stlFile().nbFacets(); ++i)
            {
                Integer aFacetId = aStl.stlFile().facetId(i);

                CGeoPolygon<Real> aPolygon;

                aPolygon.polyline().addVertex(aStl.stlFile().facet(aFacetId).vertex(0));
                aPolygon.polyline().addVertex(aStl.stlFile().facet(aFacetId).vertex(1));
                aPolygon.polyline().addVertex(aStl.stlFile().facet(aFacetId).vertex(2));
                aPolygon.close();

                aPolygon.calculateNormal();

                this->m_polygons.push_back(aPolygon);
            }
        }

        template<typename Real> 
        ENigMA::stl::CStlUtils<Real> CCsgBoolean<Real>::toStl()
        {
            ENigMA::stl::CStlUtils<Real> aStl;

            aStl.stlFile().stats.type = ENigMA::stl::FT_ASCII;

            Integer n = 0;

            for (Integer i = 0; i < static_cast<Integer> (this->m_polygons.size()); ++i)
            {
                std::vector<CGeoTriangle<Real>> sTriangles = this->m_polygons[i].triangulate();

                for (Integer j = 0; j < static_cast<Integer> (sTriangles.size()); ++j)
                {
                    ENigMA::stl::CStlFacet<Real> aFacet;

                    aFacet.addVertex(sTriangles[j].vertex(0));
                    aFacet.addVertex(sTriangles[j].vertex(1));
                    aFacet.addVertex(sTriangles[j].vertex(2));

                    aStl.addFacet(n, aFacet);
                    n++;
                }
            }

            return aStl;
        }
    }
}
