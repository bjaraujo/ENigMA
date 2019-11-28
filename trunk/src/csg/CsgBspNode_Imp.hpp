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
        CCsgBspNode<Real>::CCsgBspNode(std::vector<CGeoPolygon<Real> >& sPolygons, const Real aTolerance)
        {
            m_tolerance = aTolerance;

            m_front = NULL;
            m_back = NULL;

            if (sPolygons.size() > 0)
                this->build(sPolygons);
        }

        template<typename Real> 
        CCsgBspNode<Real>::CCsgBspNode(const Real aTolerance) 
        {
            m_tolerance = aTolerance;

            m_front = NULL;
            m_back = NULL;
        }

        template<typename Real> 
        CCsgBspNode<Real>::~CCsgBspNode() 
        {
        }

        template<typename Real> 
        void CCsgBspNode<Real>::invert()
        {
            for (Integer i = 0; i < static_cast<Integer> (this->m_polygons.size()); ++i)
                  this->m_polygons[i].invert();

            this->m_plane.invert();

            if (this->m_front)
                this->m_front->invert();

            if (this->m_back)
                this->m_back->invert();

            std::swap(this->m_front, this->m_back);
        }

        template<typename Real> 
        std::vector<CGeoPolygon<Real> > CCsgBspNode<Real>::clipPolygons(std::vector<CGeoPolygon<Real> >& sPolygons, bool removeCoplanarFront)
        {
            std::vector<CGeoPolygon<Real> > sFront;
            std::vector<CGeoPolygon<Real> > sBack;

            for (Integer i = 0; i < static_cast<Integer> (sPolygons.size()); ++i)
            {
                if (removeCoplanarFront)
                    this->splitPolygon(sPolygons[i], sBack, sBack, sFront, sBack);
                else
                    this->splitPolygon(sPolygons[i], sFront, sBack, sFront, sBack);
            }

            if (this->m_front) 
                sFront = this->m_front->clipPolygons(sFront, removeCoplanarFront);

            if (this->m_back) 
                sBack = this->m_back->clipPolygons(sBack, removeCoplanarFront);
            else
                sBack.clear();

            for (Integer i = 0; i < static_cast<Integer> (sBack.size()); ++i)
                sFront.push_back(sBack[i]);

            return sFront;
        }

        template<typename Real> 
        void CCsgBspNode<Real>::clipTo(CCsgBspNode<Real> aCsgBspNode, bool removeCoplanarFront)
        {
            this->m_polygons = aCsgBspNode.clipPolygons(this->m_polygons, removeCoplanarFront);

            if (this->m_front) 
                this->m_front->clipTo(aCsgBspNode, removeCoplanarFront);
            
            if (this->m_back) 
                this->m_back->clipTo(aCsgBspNode, removeCoplanarFront);
        }

        template<typename Real> 
        std::vector<CGeoPolygon<Real> >& CCsgBspNode<Real>::allPolygons()
        {
            if (this->m_front)
            { 
                std::vector<CGeoPolygon<Real> > allFront = this->m_front->allPolygons();

                for (Integer i = 0; i < static_cast<Integer> (allFront.size()); ++i)
                    this->m_polygons.push_back(allFront[i]);
            }

            if (this->m_back)
            {
                std::vector<CGeoPolygon<Real> > allBack = this->m_back->allPolygons();

                for (Integer i = 0; i < static_cast<Integer> (allBack.size()); ++i)
                    this->m_polygons.push_back(allBack[i]);
            }

            return this->m_polygons;
        }

        template<typename Real> 
        void CCsgBspNode<Real>::build(std::vector<CGeoPolygon<Real> >& sPolygons)
        {
            if (sPolygons.size() == 0) 
                return;

            if (sPolygons[0].polyline().nbVertices() == 0)
                return;

            Real d = sPolygons[0].normal().dot(sPolygons[0].polyline().vertex(0));
            this->m_plane.set(sPolygons[0].normal(), d);

            std::vector<CGeoPolygon<Real> > sFront;
            std::vector<CGeoPolygon<Real> > sBack;

            for (Integer i = 0; i < static_cast<Integer> (sPolygons.size()); ++i)
                this->splitPolygon(sPolygons[i], this->m_polygons, this->m_polygons, sFront, sBack);

            if (sFront.size() > 0)
            {
                if (!this->m_front)
                    this->m_front = new CCsgBspNode<Real>(m_tolerance);

                this->m_front->build(sFront);

            }

            if (sBack.size() > 0)
            {
                if (!this->m_back)
                    this->m_back = new CCsgBspNode<Real>(m_tolerance);

                this->m_back->build(sBack);

            }
        }

        template<typename Real> 
        void CCsgBspNode<Real>::splitPolygon(CGeoPolygon<Real> aPolygon, std::vector<CGeoPolygon<Real> >& sCoplanarFront, std::vector<CGeoPolygon<Real> >& sCoplanarBack, std::vector<CGeoPolygon<Real> >& sFront, std::vector<CGeoPolygon<Real> >& sBack)
        {
            const int CSG_COPLANAR = 0;
            const int CSG_FRONT = 1;
            const int CSG_BACK = 2;
            const int CSG_SPANNING = 3;

            int aPolygonType = CSG_COPLANAR;

            std::vector<int> sVertexTypes;

            for (Integer i = 0; i < aPolygon.polyline().nbVertices(); ++i)
            {
                Real t = this->m_plane.normal().dot(aPolygon.polyline().vertex(i)) - this->m_plane.d();

                int aVertexType = (t < -m_tolerance) ? CSG_BACK : (t > m_tolerance) ? CSG_FRONT : CSG_COPLANAR;

                aPolygonType |= aVertexType;
                sVertexTypes.push_back(aVertexType);

            }

            switch (aPolygonType)
            {
            case CSG_COPLANAR:

                (this->m_plane.normal().dot(aPolygon.normal()) > 0 ? sCoplanarFront : sCoplanarBack).push_back(aPolygon);
                break;

            case CSG_FRONT:

                sFront.push_back(aPolygon);
                break;

            case CSG_BACK:

                sBack.push_back(aPolygon);
                break;

            case CSG_SPANNING:

                CGeoPolygon<Real> aPolygonFront;
                CGeoPolygon<Real> aPolygonBack;

                for (Integer i = 0; i < aPolygon.polyline().nbVertices(); ++i)
                {
                    Integer j = (i + 1) % aPolygon.polyline().nbVertices();

                    int ti = sVertexTypes[i];
                    int tj = sVertexTypes[j];

                    CGeoCoordinate<Real> vi = aPolygon.polyline().vertex(i);
                    CGeoCoordinate<Real> vj = aPolygon.polyline().vertex(j);

                    if (ti != CSG_BACK)
                        aPolygonFront.polyline().addVertex(vi);

                    if (ti != CSG_FRONT)
                        aPolygonBack.polyline().addVertex(vi);

                    if ((ti | tj) == CSG_SPANNING)
                    {
                        Real t = (this->m_plane.d() - this->m_plane.normal().dot(vi)) / this->m_plane.normal().dot(vj - vi);
                        CGeoCoordinate<Real> v = vi + (vj - vi) * t;

                        aPolygonFront.polyline().addVertex(v);
                        aPolygonBack.polyline().addVertex(v);
                    }

                }

                if (aPolygonFront.polyline().nbVertices() > 2)
                {
                    aPolygonFront.polyline().close();
                    aPolygonFront.normal() = aPolygon.normal();

                    sFront.push_back(aPolygonFront);
                }

                if (aPolygonBack.polyline().nbVertices() > 2)
                {
                    aPolygonBack.polyline().close();
                    aPolygonBack.normal() = aPolygon.normal();

                    sBack.push_back(aPolygonBack);
                }

                break;

            }
        }

        template<typename Real> 
        void CCsgBspNode<Real>::clear()
        {
            m_polygons.clear();

            if (this->m_back)
            {
                this->m_back->clear();
                delete this->m_back;
            }

            if (this->m_front)
            {
                this->m_front->clear();
                delete this->m_front;
            }
        }
    }
}
