// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#include "gtest/gtest.h"

#include "TypeDef.hpp"

#include "GeoPolygon.hpp"

using namespace ENigMA::geometry;

class CTestGeoPolygon : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoPolygon, area1) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    aPolygon.calculateNormal();

    EXPECT_EQ(1.0, aPolygon.normal().z());

    aPolygon.calculateArea();

    EXPECT_EQ(1.0, aPolygon.area());

}

TEST_F(CTestGeoPolygon, area2) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(3.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(3.0, 3.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 3.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    aPolygon.calculateNormal();

    EXPECT_EQ(1.0, aPolygon.normal().z());

    aPolygon.calculateArea();

    EXPECT_EQ(9.0, aPolygon.area());

}

TEST_F(CTestGeoPolygon, centroid) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    aPolygon.calculateCentroid();

    EXPECT_EQ(0.5, aPolygon.centroid().x());
    EXPECT_EQ(0.5, aPolygon.centroid().y());
    EXPECT_EQ(0.0, aPolygon.centroid().z());

}

TEST_F(CTestGeoPolygon, clip1) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(2.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 0.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, 1.3);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane, 1E-6);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(1.3, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip2) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 0.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, 1.0);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(1.0, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip3) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 1.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, sqrt(2.0) / 2.0);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(0.5, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip4) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 1.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, sqrt(2.0));

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(1.0, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip5) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(2.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex6(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex6);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 0.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, 1.0);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(1.0, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip6) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(2.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(2.0, 2.0, 0.0);
    CGeoCoordinate<decimal> aVertex6(0.0, 2.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex6);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 1.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, 3.0 * sqrt(2.0) / 2.0);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(2.5, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip7) {

    CGeoCoordinate<decimal> aVertex1(+0.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex3(+2.0, +2.0, +0.0);
    CGeoCoordinate<decimal> aVertex4(+1.0, +3.0, +0.0);
    CGeoCoordinate<decimal> aVertex5(+0.0, +3.0, +0.0);
    CGeoCoordinate<decimal> aVertex6(-0.5, +2.5, +0.0);
    CGeoCoordinate<decimal> aVertex7(-0.5, +1.5, +0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex6);
    aPolyline.addVertex(aVertex7);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 0.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, 1.5);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(3.5, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip8) {

    CGeoCoordinate<decimal> aVertex1(+0.0, +1.0, +1.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +1.0, +1.0);
    CGeoCoordinate<decimal> aVertex3(+1.5, +1.5, +1.0);
    CGeoCoordinate<decimal> aVertex4(+1.5, +2.5, +1.0);
    CGeoCoordinate<decimal> aVertex5(+1.0, +3.0, +1.0);
    CGeoCoordinate<decimal> aVertex6(+0.0, +3.0, +1.0);
    CGeoCoordinate<decimal> aVertex7(-0.5, +2.5, +1.0);
    CGeoCoordinate<decimal> aVertex8(-0.5, +1.5, +1.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex6);
    aPolyline.addVertex(aVertex7);
    aPolyline.addVertex(aVertex8);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 1.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, 5.0 * sqrt(2.0) / 4.0);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(3.5 / 2.0, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip9) {

    CGeoCoordinate<decimal> aVertex1(+0.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex3(+1.0, +0.2, +0.0);
    CGeoCoordinate<decimal> aVertex4(+0.2, +0.2, +0.0);
    CGeoCoordinate<decimal> aVertex5(+0.2, +0.8, +0.0);
    CGeoCoordinate<decimal> aVertex6(+1.0, +0.8, +0.0);
    CGeoCoordinate<decimal> aVertex7(+1.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex8(+0.0, +1.0, +0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex5);
    aPolyline.addVertex(aVertex6);
    aPolyline.addVertex(aVertex7);
    aPolyline.addVertex(aVertex8);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    CGeoNormal<decimal> aNormal(1.0, 0.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, 0.4);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(0.28, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, clip10) {

    CGeoCoordinate<decimal> aVertex1(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex3(+1.0, +1.0, +1.0);
    CGeoCoordinate<decimal> aVertex4(+1.0, +0.0, +1.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex4);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex1);
    aPolyline.close();

    CGeoPolygon<decimal> aPolygon;

    aPolygon.setPolyline(aPolyline);

    aPolygon.calculateArea();

    EXPECT_NEAR(1.0, aPolygon.area(), 1E-6);

    CGeoNormal<decimal> aNormal(1.0, 1.0, 0.0);

    CGeoPlane<decimal> aPlane(aNormal, sqrt(2.0) / 2.0);

    CGeoPolygon<decimal> aNewPolygon = aPolygon.clip(aPlane);

    aNewPolygon.calculateArea();

    EXPECT_NEAR(0.0, aNewPolygon.area(), 1E-6);

}

TEST_F(CTestGeoPolygon, triangulate) {

    CGeoCoordinate<decimal> aVertex1(+0.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex3(+2.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex4(+2.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex5(+1.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex6(+0.0, +1.0, +0.0);

    CGeoPolygon<decimal> aPolygon;

    aPolygon.polyline().addVertex(aVertex1);
    aPolygon.polyline().addVertex(aVertex2);
    aPolygon.polyline().addVertex(aVertex3);
    aPolygon.polyline().addVertex(aVertex4);
    aPolygon.polyline().addVertex(aVertex5);
    aPolygon.polyline().addVertex(aVertex6);
    aPolygon.polyline().close();

    std::vector<CGeoTriangle<decimal> > sTriangles = aPolygon.triangulate();

    EXPECT_EQ(4, sTriangles.size());

}

TEST_F(CTestGeoPolygon, removeCollinear) {

    CGeoCoordinate<decimal> aVertex1(+0.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex3(+2.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex4(+2.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex5(+1.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex6(+0.0, +1.0, +0.0);

    CGeoPolygon<decimal> aPolygon;

    aPolygon.polyline().addVertex(aVertex1);
    aPolygon.polyline().addVertex(aVertex2);
    aPolygon.polyline().addVertex(aVertex3);
    aPolygon.polyline().addVertex(aVertex4);
    aPolygon.polyline().addVertex(aVertex5);
    aPolygon.polyline().addVertex(aVertex6);
    aPolygon.polyline().close();

    aPolygon.polyline().removeCollinear(1E-6);

    EXPECT_EQ(5, aPolygon.polyline().nbVertices());

}
