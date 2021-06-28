// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iomanip>
#include <iostream>

#include <time.h>

#include "GeoPolyhedron.hpp"

using namespace ENigMA::geometry;

int main(int argc, char* argv[])
{
    CGeoCoordinate<float> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<float> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<float> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<float> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<float> aVertex5(0.0, 0.0, -1.0);
    CGeoCoordinate<float> aVertex6(1.0, 0.0, -1.0);
    CGeoCoordinate<float> aVertex7(1.0, 1.0, -1.0);
    CGeoCoordinate<float> aVertex8(0.0, 1.0, -1.0);

    CGeoHexahedron<float> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CGeoPolyhedron<float> aPolyhedron(aHexahedron);

    CGeoNormal<float> aNormal(0.9407, 0.2822, 0.1881);
    CGeoPlane<float> aPlane(aNormal, 0.0);

    CGeoPolygon<float> aNewPolygon;
    Integer aNewPolygonId = 999;

    float volumeFractionAct;
    Integer nIterations;

    CGeoPolyhedron<float> aNewPolyhedron = aPolyhedron.clip(aNewPolygon, aNewPolygonId, aPlane, 0.4, volumeFractionAct, nIterations, 50, 1E-9, 1E-6);

    std::cout << "Solution: " << aPlane.d() << std::endl;

    return 0;
}
