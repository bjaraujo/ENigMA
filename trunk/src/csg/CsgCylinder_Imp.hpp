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
        CGeoCoordinate<Real> CCsgCylinder<Real>::point(CGeoCoordinate<Real>& aStart, CGeoVector<Real>& axisX, CGeoVector<Real>& axisY, CGeoVector<Real>& axisZ, CGeoVector<Real>& aRay, Real aRadius, Real stack, Real slice, Real normalBlend)
        {
            
            const Real pi = std::acos(-1.0);

            Real angle = slice * pi * 2;

            CGeoVector<Real> out = axisX * cos(angle) + axisY * sin(angle);
        
            return aStart + aRay * stack + out * aRadius;
        }

        template<typename Real> 
        ENigMA::csg::CCsgBoolean<Real> CCsgCylinder<Real>::create(CGeoCoordinate<Real>& aStart, CGeoCoordinate<Real>& anEnd, Real aRadius, Integer nSlices)
        {
            std::vector<CGeoPolygon<Real> > sPolygons;

            CGeoVector<Real> aRay = anEnd - aStart;
            CGeoVector<Real> axisZ = aRay;
            axisZ.normalize();

            CGeoVector<Real> v;

            if (fabs(axisZ.y()) > 0.5)
                v << 1, 0, 0;
            else
                v << 0, 1, 0;

            CGeoVector<Real> axisX = v.cross(axisZ);
            axisX.normalize();
            
            CGeoVector<Real> axisY = axisX.cross(axisZ);
            axisY.normalize();

            for (Integer i = 0; i < nSlices; ++i) 
            {
                
                Real t0 = ((Real) i) / nSlices; 
                Real t1 = ((Real) (i + 1)) / nSlices;

                CGeoPolygon<Real> aPolygon1;

                aPolygon1.polyline().addVertex(aStart);

                CGeoCoordinate<Real> aVertex1 = CCsgCylinder<Real>::point(aStart, axisX, axisY, axisZ, aRay, aRadius, 0, t0, -1);
                CGeoCoordinate<Real> aVertex2 = CCsgCylinder<Real>::point(aStart, axisX, axisY, axisZ, aRay, aRadius, 0, t1, -1);

                aPolygon1.polyline().addVertex(aVertex1);
                aPolygon1.polyline().addVertex(aVertex2);

                aPolygon1.close();
                aPolygon1.calculateNormal();

                sPolygons.push_back(aPolygon1);

                CGeoPolygon<Real> aPolygon2;

                CGeoCoordinate<Real> aVertex3 = CCsgCylinder<Real>::point(aStart, axisX, axisY, axisZ, aRay, aRadius, 0, t1, 0);
                CGeoCoordinate<Real> aVertex4 = CCsgCylinder<Real>::point(aStart, axisX, axisY, axisZ, aRay, aRadius, 0, t0, 0);
                CGeoCoordinate<Real> aVertex5 = CCsgCylinder<Real>::point(aStart, axisX, axisY, axisZ, aRay, aRadius, 1, t0, 0);
                CGeoCoordinate<Real> aVertex6 = CCsgCylinder<Real>::point(aStart, axisX, axisY, axisZ, aRay, aRadius, 1, t1, 0);

                aPolygon2.polyline().addVertex(aVertex3);
                aPolygon2.polyline().addVertex(aVertex4);
                aPolygon2.polyline().addVertex(aVertex5);
                aPolygon2.polyline().addVertex(aVertex6);

                aPolygon2.close();
                aPolygon2.calculateNormal();

                sPolygons.push_back(aPolygon2);

                CGeoPolygon<Real> aPolygon3;

                aPolygon3.polyline().addVertex(anEnd);

                CGeoCoordinate<Real> aVertex7 = CCsgCylinder<Real>::point(aStart, axisX, axisY, axisZ, aRay, aRadius, 1, t1, 1);
                CGeoCoordinate<Real> aVertex8 = CCsgCylinder<Real>::point(aStart, axisX, axisY, axisZ, aRay, aRadius, 1, t0, 1);

                aPolygon3.polyline().addVertex(aVertex7);
                aPolygon3.polyline().addVertex(aVertex8);

                aPolygon3.close();
                aPolygon3.calculateNormal();

                sPolygons.push_back(aPolygon3);

            }

            ENigMA::csg::CCsgBoolean<Real> result(sPolygons);
            return result;
        }
    }
}
