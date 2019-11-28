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
        ENigMA::csg::CCsgBoolean<Real> CCsgCube<Real>::create(CGeoCoordinate<Real> aCenter, CGeoVector<Real> aRadius)
        {
            std::vector<CGeoCoordinate<Real>> sVertices;

            for (Integer i = 0; i < 8; ++i)
            {
                Real x = aCenter.x() + aRadius.x() * (2 * !!(i & 1) - 1);
                Real y = aCenter.y() + aRadius.y() * (2 * !!(i & 2) - 1);
                Real z = aCenter.z() + aRadius.z() * (2 * !!(i & 4) - 1);

                CGeoCoordinate<Real> aVertex(x, y, z);

                sVertices.push_back(aVertex);

            }

            std::vector<CGeoPolygon<Real>> sPolygons;

            // face 1
            CGeoPolygon<Real> aPolygon1;
            aPolygon1.polyline().addVertex(sVertices[0]);
            aPolygon1.polyline().addVertex(sVertices[4]);
            aPolygon1.polyline().addVertex(sVertices[6]);
            aPolygon1.polyline().addVertex(sVertices[2]);
            aPolygon1.polyline().close();
            aPolygon1.calculateNormal();
            sPolygons.push_back(aPolygon1);

            // face 2
            CGeoPolygon<Real> aPolygon2;
            aPolygon2.polyline().addVertex(sVertices[1]);
            aPolygon2.polyline().addVertex(sVertices[3]);
            aPolygon2.polyline().addVertex(sVertices[7]);
            aPolygon2.polyline().addVertex(sVertices[5]);
            aPolygon2.polyline().close();
            aPolygon2.calculateNormal();
            sPolygons.push_back(aPolygon2);

            // face 3
            CGeoPolygon<Real> aPolygon3;
            aPolygon3.polyline().addVertex(sVertices[0]);
            aPolygon3.polyline().addVertex(sVertices[1]);
            aPolygon3.polyline().addVertex(sVertices[5]);
            aPolygon3.polyline().addVertex(sVertices[4]);
            aPolygon3.polyline().close();
            aPolygon3.calculateNormal();
            sPolygons.push_back(aPolygon3);

            // face 4
            CGeoPolygon<Real> aPolygon4;
            aPolygon4.polyline().addVertex(sVertices[2]);
            aPolygon4.polyline().addVertex(sVertices[6]);
            aPolygon4.polyline().addVertex(sVertices[7]);
            aPolygon4.polyline().addVertex(sVertices[3]);
            aPolygon4.polyline().close();
            aPolygon4.calculateNormal();
            sPolygons.push_back(aPolygon4);

            // face 5
            CGeoPolygon<Real> aPolygon5;
            aPolygon5.polyline().addVertex(sVertices[0]);
            aPolygon5.polyline().addVertex(sVertices[2]);
            aPolygon5.polyline().addVertex(sVertices[3]);
            aPolygon5.polyline().addVertex(sVertices[1]);
            aPolygon5.polyline().close();
            aPolygon5.calculateNormal();
            sPolygons.push_back(aPolygon5);

            // face 6
            CGeoPolygon<Real> aPolygon6;
            aPolygon6.polyline().addVertex(sVertices[4]);
            aPolygon6.polyline().addVertex(sVertices[5]);
            aPolygon6.polyline().addVertex(sVertices[7]);
            aPolygon6.polyline().addVertex(sVertices[6]);
            aPolygon6.polyline().close();
            aPolygon6.calculateNormal();
            sPolygons.push_back(aPolygon6);

            ENigMA::csg::CCsgBoolean<Real> result(sPolygons);
            return result;
        }
    }
}
