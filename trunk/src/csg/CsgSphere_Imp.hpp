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
        template <typename Real>
        CGeoCoordinate<Real> CCsgSphere<Real>::vertex(CGeoCoordinate<Real>& aCenter, Real aRadius, Real theta, Real phi)
        {
            const Real pi = std::acos(-1.0);

            theta *= pi * 2;
            phi *= pi;

            CGeoVector<Real> dir(cos(theta) * sin(phi), cos(phi), sin(theta) * sin(phi));

            return aCenter + (dir * aRadius);
        }

        template <typename Real>
        ENigMA::csg::CCsgBoolean<Real> CCsgSphere<Real>::create(CGeoCoordinate<Real>& aCenter, Real aRadius, Integer nSlices, Integer nStacks)
        {
            std::vector<CGeoPolygon<Real>> sPolygons;

            for (Integer i = 0; i < nSlices; ++i)
            {
                for (Integer j = 0; j < nStacks; ++j)
                {
                    CGeoPolygon<Real> aPolygon;

                    CGeoCoordinate<Real> aVertex1 = CCsgSphere<Real>::vertex(aCenter, aRadius, ((Real)i) / nSlices, ((Real)j) / nStacks);
                    CGeoCoordinate<Real> aVertex2 = CCsgSphere<Real>::vertex(aCenter, aRadius, ((Real)(i + 1)) / nSlices, ((Real)j) / nStacks);
                    CGeoCoordinate<Real> aVertex3 = CCsgSphere<Real>::vertex(aCenter, aRadius, ((Real)(i + 1)) / nSlices, ((Real)(j + 1)) / nStacks);
                    CGeoCoordinate<Real> aVertex4 = CCsgSphere<Real>::vertex(aCenter, aRadius, ((Real)i) / nSlices, ((Real)(j + 1)) / nStacks);

                    aPolygon.polyline().addVertex(aVertex1);

                    if (j > 0)
                        aPolygon.polyline().addVertex(aVertex2);
                    if (j < nStacks - 1)
                        aPolygon.polyline().addVertex(aVertex3);

                    aPolygon.polyline().addVertex(aVertex4);

                    aPolygon.polyline().close();
                    aPolygon.calculateNormal();

                    sPolygons.push_back(aPolygon);
                }
            }

            ENigMA::csg::CCsgBoolean<Real> result(sPolygons);
            return result;
        }
    }
}
