// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "gtest/gtest.h"

#include "TypeDef.hpp"

#include "GeoPolyhedron.hpp"

using namespace ENigMA::geometry;

class CTestGeoPolyhedron : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoPolyhedron, centroid) {

    CGeoPolyline<Decimal> aPolyline;
    CGeoPolygon<Decimal> aPolygon;
    CGeoPolyhedron<Decimal> aPolyhedron;

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, 1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, 1.0);
    CGeoCoordinate<Decimal> aVertex7(1.0, 1.0, 1.0);
    CGeoCoordinate<Decimal> aVertex8(0.0, 1.0, 1.0);

    // face 1
    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex2);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(0, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 2
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex6);
    aPolyline.addVertex(aVertex7);
    aPolyline.addVertex(aVertex8);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(1, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 3
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex6);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(2, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 4
    aPolyline.addVertex(aVertex8);
    aPolyline.addVertex(aVertex7);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(3, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 5
    aPolyline.addVertex(aVertex6);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex7);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(4, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 6
    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex8);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(5, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    aPolyhedron.calculateCentroid(true);

    EXPECT_EQ(0.5, aPolyhedron.centroid().x());
    EXPECT_EQ(0.5, aPolyhedron.centroid().y());
    EXPECT_EQ(0.5, aPolyhedron.centroid().z());

}

TEST_F(CTestGeoPolyhedron, volume1) {

    CGeoPolyline<Decimal> aPolyline;
    CGeoPolygon<Decimal> aPolygon;
    CGeoPolyhedron<Decimal> aPolyhedron;

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex7(1.0, 1.0, -1.0);
    CGeoCoordinate<Decimal> aVertex8(0.0, 1.0, -1.0);

    // face 1
    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(0, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 2
    aPolyline.addVertex(aVertex8);
    aPolyline.addVertex(aVertex7);
    aPolyline.addVertex(aVertex6);
    aPolyline.addVertex(aVertex5);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(1, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 3
    aPolyline.addVertex(aVertex6);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex5);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(2, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 4
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex7);
    aPolyline.addVertex(aVertex8);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(3, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 5
    aPolyline.addVertex(aVertex7);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex6);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(4, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 6
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex8);
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex1);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);
    aPolyhedron.addPolygon(5, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    aPolyhedron.calculateSurfaceArea(true);

    EXPECT_EQ(6.0, aPolyhedron.surfaceArea());

    aPolyhedron.calculateVolume(true);

    EXPECT_EQ(1.0, aPolyhedron.volume());

}

TEST_F(CTestGeoPolyhedron, volume2) {

    Decimal px, py, pz;

    px = +2.0;
    py = +5.0;
    pz = -20.0;

    CGeoPolyline<Decimal> aPolyline;
    CGeoPolygon<Decimal> aPolygon;
    CGeoPolyhedron<Decimal> aPolyhedron;

    CGeoCoordinate<Decimal> aVertex1(px + 0.0, py + 0.0, pz + 0.0);
    CGeoCoordinate<Decimal> aVertex2(px + 1.0, py + 0.0, pz + 0.0);
    CGeoCoordinate<Decimal> aVertex3(px + 0.0, py + 1.0, pz + 0.0);
    CGeoCoordinate<Decimal> aVertex4(px + 0.0, py + 0.0, pz - 1.0);

    // face 1
    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);

    aPolygon.calculateNormal(true);
    EXPECT_NEAR(1.0, aPolygon.normal().z(), 1E-6);

    aPolygon.calculateArea(true);
    EXPECT_NEAR(0.5, aPolygon.area(), 1E-6);

    aPolyhedron.addPolygon(0, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 2
    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex2);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);

    aPolygon.calculateNormal(true);
    EXPECT_NEAR(-1.0, aPolygon.normal().y(), 1E-6);

    aPolygon.calculateArea(true);
    EXPECT_NEAR(0.5, aPolygon.area(), 1E-6);

    aPolyhedron.addPolygon(1, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 3
    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);

    aPolygon.calculateNormal(true);
    EXPECT_NEAR(-1.0, aPolygon.normal().x(), 1E-6);

    aPolygon.calculateArea(true);
    EXPECT_NEAR(0.5, aPolygon.area(), 1E-6);

    aPolyhedron.addPolygon(2, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    // face 4
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex3);
    aPolyline.close();

    aPolygon.setPolyline(aPolyline);

    aPolygon.calculateArea(true);
    EXPECT_NEAR(sqrt(1.5) * sqrt(2.0) / 2.0, aPolygon.area(), 1E-6);

    aPolyhedron.addPolygon(3, aPolygon);
    aPolyline.reset();
    aPolygon.reset();

    aPolyhedron.calculateSurfaceArea(true);

    EXPECT_NEAR(sqrt(1.5) * sqrt(2.0) / 2.0 + 1.5, aPolyhedron.surfaceArea(), 1E-6);

    aPolyhedron.calculateVolume(true);

    EXPECT_NEAR(1.0/6.0, aPolyhedron.volume(), 1E-6);

}

TEST_F(CTestGeoPolyhedron, volume3) {

    Decimal px, py, pz;

    px = 2.0;
    py = 5.0;
    pz = 20.0;

    CGeoCoordinate<Decimal> aVertex1(px + 0.0, py + 0.0, pz + 0.0);
    CGeoCoordinate<Decimal> aVertex2(px + 1.0, py + 0.0, pz + 0.0);
    CGeoCoordinate<Decimal> aVertex3(px + 1.0, py + 1.0, pz + 0.0);
    CGeoCoordinate<Decimal> aVertex4(px + 0.0, py + 1.0, pz + 0.0);
    CGeoCoordinate<Decimal> aVertex5(px + 0.0, py + 0.0, pz - 1.0);
    CGeoCoordinate<Decimal> aVertex6(px + 1.0, py + 0.0, pz - 1.0);
    CGeoCoordinate<Decimal> aVertex7(px + 1.0, py + 1.0, pz - 1.0);
    CGeoCoordinate<Decimal> aVertex8(px + 0.0, py + 1.0, pz - 1.0);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    aHexahedron.calculateVolume();

    CGeoPolyhedron<Decimal> aPolyhedron(aHexahedron);

    aPolyhedron.calculateSurfaceArea(true);

    EXPECT_EQ(6.0, aPolyhedron.surfaceArea());

    aPolyhedron.calculateVolume(true);

    EXPECT_EQ(1.0, aPolyhedron.volume());

    EXPECT_NEAR(aPolyhedron.volume(), aHexahedron.volume(), 1E-6);

}

TEST_F(CTestGeoPolyhedron, volume4) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 0.0, -1.0);

    CGeoTetrahedron<Decimal> aTetrahedron;

    aTetrahedron.addVertex(aVertex1);
    aTetrahedron.addVertex(aVertex2);
    aTetrahedron.addVertex(aVertex3);
    aTetrahedron.addVertex(aVertex4);

    aTetrahedron.calculateVolume();

    EXPECT_NEAR(1.0/6.0, aTetrahedron.volume(), 1E-12);

    CGeoPolyhedron<Decimal> aPolyhedron(aTetrahedron);

    aPolyhedron.calculateVolume();

    EXPECT_NEAR(1.0/6.0, aPolyhedron.volume(), 1E-6);

}

TEST_F(CTestGeoPolyhedron, volume5) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex5(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 1.0, -1.0);

    CGeoTriangularPrism<Decimal> aTriangularPrism;

    aTriangularPrism.addVertex(aVertex1);
    aTriangularPrism.addVertex(aVertex2);
    aTriangularPrism.addVertex(aVertex3);
    aTriangularPrism.addVertex(aVertex4);
    aTriangularPrism.addVertex(aVertex5);
    aTriangularPrism.addVertex(aVertex6);

    aTriangularPrism.calculateVolume();

    EXPECT_EQ(0.5, aTriangularPrism.volume());

    CGeoPolyhedron<Decimal> aPolyhedron(aTriangularPrism);

    aPolyhedron.calculateVolume();

    EXPECT_NEAR(0.5, aPolyhedron.volume(), 1E-12);

}

TEST_F(CTestGeoPolyhedron, volume6) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex7(1.0, 1.0, -1.0);
    CGeoCoordinate<Decimal> aVertex8(0.0, 1.0, -1.0);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    aHexahedron.calculateVolume();

    EXPECT_NEAR(1.0, aHexahedron.volume(), 1E-12);

    CGeoPolyhedron<Decimal> aPolyhedron(aHexahedron);

    aPolyhedron.calculateVolume();

    EXPECT_NEAR(1.0, aPolyhedron.volume(), 1E-12);

}

TEST_F(CTestGeoPolyhedron, clip1) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex7(1.0, 1.0, -1.0);
    CGeoCoordinate<Decimal> aVertex8(0.0, 1.0, -1.0);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    aHexahedron.calculateVolume();

    CGeoPolyhedron<Decimal> aPolyhedron(aHexahedron);

    aPolyhedron.calculateSurfaceArea(true);

    EXPECT_EQ(6.0, aPolyhedron.surfaceArea());

    aPolyhedron.calculateVolume(true);

    EXPECT_EQ(1.0, aPolyhedron.volume());

    CGeoNormal<Decimal> aNormal(1.0, 0.0, 0.0);

    CGeoPlane<Decimal> aPlane(aNormal, 0.5);

    CGeoPolygon<Decimal> aNewPolygon;
    Integer aNewPolygonId = 999;

    CGeoPolyhedron<Decimal> aNewPolyhedron = aPolyhedron.clip(aNewPolygon, aNewPolygonId, aPlane, 1E-6);

    EXPECT_EQ(4, aNewPolygon.polyline().nbLines());

    EXPECT_EQ(6, aNewPolyhedron.nbPolygons());

    aNewPolyhedron.calculateSurfaceArea(true);

    EXPECT_EQ(4.0, aNewPolyhedron.surfaceArea());

    aNewPolyhedron.calculateVolume(true);

    EXPECT_NEAR(0.5, aNewPolyhedron.volume(), 1E-6);

    EXPECT_NEAR(aPolyhedron.volume(), aHexahedron.volume(), 1E-6);

}

TEST_F(CTestGeoPolyhedron, clip2) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex7(1.0, 1.0, -1.0);
    CGeoCoordinate<Decimal> aVertex8(0.0, 1.0, -1.0);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CGeoPolyhedron<Decimal> aPolyhedron(aHexahedron);

    EXPECT_EQ(6, aPolyhedron.nbPolygons());

    aPolyhedron.calculateSurfaceArea(true);

    EXPECT_EQ(6.0, aPolyhedron.surfaceArea());

    aPolyhedron.calculateVolume(true);

    EXPECT_EQ(1.0, aPolyhedron.volume());

    CGeoNormal<Decimal> aNormal(1.0, 1.0, 0.0);

    CGeoPlane<Decimal> aPlane(aNormal, sqrt(2.0) / 2.0);

    CGeoPolygon<Decimal> aNewPolygon;
    Integer aNewPolygonId = 999;

    CGeoPolyhedron<Decimal> aNewPolyhedron = aPolyhedron.clip(aNewPolygon, aNewPolygonId, aPlane, 1E-6);

    EXPECT_EQ(5, aNewPolyhedron.nbPolygons());

    aNewPolyhedron.calculateSurfaceArea(true);

    EXPECT_NEAR(3.0 + sqrt(2.0), aNewPolyhedron.surfaceArea(), 1E-6);

    aNewPolyhedron.calculateVolume(true);

    EXPECT_NEAR(0.5, aNewPolyhedron.volume(), 1E-6);

}

TEST_F(CTestGeoPolyhedron, clip3) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex5(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 1.0, -1.0);

    CGeoTriangularPrism<Decimal> aTriangularPrism;

    aTriangularPrism.addVertex(aVertex1);
    aTriangularPrism.addVertex(aVertex2);
    aTriangularPrism.addVertex(aVertex3);
    aTriangularPrism.addVertex(aVertex4);
    aTriangularPrism.addVertex(aVertex5);
    aTriangularPrism.addVertex(aVertex6);

    aTriangularPrism.calculateVolume();

    EXPECT_EQ(0.5, aTriangularPrism.volume());

    CGeoPolyhedron<Decimal> aPolyhedron(aTriangularPrism);

    aPolyhedron.calculateVolume();

    EXPECT_NEAR(0.5, aPolyhedron.volume(), 1E-12);

    CGeoNormal<Decimal> aNormal(1.0, 1.0, 0.0);

    CGeoPlane<Decimal> aPlane(aNormal, sqrt(2.0) / 2.0);

    CGeoPolygon<Decimal> aNewPolygon;
    Integer aNewPolygonId = 999;

    CGeoPolyhedron<Decimal> aNewPolyhedron = aPolyhedron.clip(aNewPolygon, aNewPolygonId, aPlane, 1E-6);

    aNewPolyhedron.calculateVolume(true);

    EXPECT_NEAR(0.25, aNewPolyhedron.volume(), 1E-6);

}

TEST_F(CTestGeoPolyhedron, clip4) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex7(1.0, 1.0, -1.0);
    CGeoCoordinate<Decimal> aVertex8(0.0, 1.0, -1.0);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CGeoPolyhedron<Decimal> aPolyhedron(aHexahedron);

    EXPECT_EQ(6, aPolyhedron.nbPolygons());

    CGeoNormal<Decimal> aNormal(0.9407, 0.2822, 0.1881);
    CGeoPlane<Decimal> aPlane(aNormal, 0.0);

    CGeoPolygon<Decimal> aNewPolygon;
    Integer aNewPolygonId = 999;

    Decimal volumeFractionAct;
    Integer nIterations;

    CGeoPolyhedron<Decimal> aNewPolyhedron = aPolyhedron.clip(aNewPolygon, aNewPolygonId, aPlane, 0.4, volumeFractionAct, nIterations, 50, 1E-9, 1E-6);

    EXPECT_EQ(5, nIterations);

    EXPECT_NEAR(0.4233437716, aPlane.d(), 1E-6);

    aNewPolyhedron.calculateVolume(true);

    EXPECT_NEAR(0.4, aNewPolyhedron.volume(), 1E-6);

}
