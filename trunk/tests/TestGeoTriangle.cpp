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

#include "GeoTriangle.hpp"

using namespace ENigMA::geometry;

class CTestGeoTriangle : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoTriangle, area) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);

    CGeoTriangle<decimal> aTriangle;

    aTriangle.addVertex(aVertex1);
    aTriangle.addVertex(aVertex2);
    aTriangle.addVertex(aVertex3);

    aTriangle.calculateArea();

    EXPECT_EQ(0.5, aTriangle.area());

}

TEST_F(CTestGeoTriangle, contains) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.5, 1.0, 0.0);

    CGeoTriangle<decimal> aTriangle;

    aTriangle.addVertex(aVertex1);
    aTriangle.addVertex(aVertex2);
    aTriangle.addVertex(aVertex3);

    CGeoCoordinate<decimal> aPoint(0.5, 0.25, 0.0);

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aTriangle.contains(aPoint, anIntersectionType));

}

TEST_F(CTestGeoTriangle, intersect1) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(0.0, 0.125, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 0.0, 0.125);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.125, 0.0);
    CGeoCoordinate<decimal> aVertex6(0.0, 0.0, 0.125);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-6);

    EXPECT_TRUE(anIntersectionType == IT_COINCIDENT);

}

TEST_F(CTestGeoTriangle, intersect2) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(0.0, 0.125, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 0.0, 0.125);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.125, 0.0);
    CGeoCoordinate<decimal> aVertex6(0.0, 0.0, -0.125);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-6);

    EXPECT_TRUE(anIntersectionType == IT_EDGE);

}

TEST_F(CTestGeoTriangle, intersect3) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(0.0, 0.125, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 0.0, 0.125);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.125, -0.125);
    CGeoCoordinate<decimal> aVertex6(0.0, 0.0, -0.125);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-6);

    EXPECT_TRUE(anIntersectionType == IT_VERTEX);

}

TEST_F(CTestGeoTriangle, intersect4) {

    CGeoCoordinate<decimal> aVertex1(0, 0.25, 0);
    CGeoCoordinate<decimal> aVertex2(0.11164981918867756, 0.13740739494625337, 0.13740739494625337);
    CGeoCoordinate<decimal> aVertex3(0.125, 0.125, 0);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(0, 0.25, 0);
    CGeoCoordinate<decimal> aVertex5(0.125, 0.125, 0);
    CGeoCoordinate<decimal> aVertex6(0.11164981918867756, 0.13740739494625337, 0.13740739494625337);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-7);

    EXPECT_TRUE(anIntersectionType == IT_COINCIDENT);

}

TEST_F(CTestGeoTriangle, intersect5) {

    CGeoCoordinate<decimal> aVertex1(0.125, 0.25, 0.25);
    CGeoCoordinate<decimal> aVertex2(0, 0.125, 0.125);
    CGeoCoordinate<decimal> aVertex3(0.083333333333333329, 0.13096452442887269, 0.12500000000000000);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(0.125, 0.25, 0.25);
    CGeoCoordinate<decimal> aVertex5(0.125, 0.25, 0);
    CGeoCoordinate<decimal> aVertex6(0.083333333333333329, 0.13096452442887269, 0.12500000000000000);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-7);

    EXPECT_TRUE(anIntersectionType == IT_EDGE);

}

TEST_F(CTestGeoTriangle, intersect6) {

    CGeoCoordinate<decimal> aVertex1(0, 0, 0);
    CGeoCoordinate<decimal> aVertex2(1, 0, 0);
    CGeoCoordinate<decimal> aVertex3(1, 1, 0);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(0, 0, 0);
    CGeoCoordinate<decimal> aVertex5(0, 1, 0);
    CGeoCoordinate<decimal> aVertex6(1, 0, 0);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-7);

    EXPECT_TRUE(anIntersectionType == IT_SWAP);

}

TEST_F(CTestGeoTriangle, intersect7) {

    CGeoCoordinate<decimal> aVertex1(-6.672, +0.503, +2.394);
    CGeoCoordinate<decimal> aVertex2(-5.880, +1.454, +2.296);
    CGeoCoordinate<decimal> aVertex3(-5.880, +0.190, +1.099);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(-6.672, +0.503, +2.394);
    CGeoCoordinate<decimal> aVertex5(-6.672, +1.156, +1.756);
    CGeoCoordinate<decimal> aVertex6(-5.419, +1.037, +1.945);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-3);

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

    aTriangle2.intersects(aTriangle1, anIntersectionType, 1E-3);

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

}

TEST_F(CTestGeoTriangle, intersect8) {

    CGeoCoordinate<decimal> aVertex1(-6.100E-002, +6.667E-001, +2.233E-002);
    CGeoCoordinate<decimal> aVertex2(+0.000E+000, +7.500E-001, +2.500E-001);
    CGeoCoordinate<decimal> aVertex3(+2.500E-001, +7.500E-001, +0.000E+000);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(+0.000E+000, +5.000E-001, +0.000E+000);
    CGeoCoordinate<decimal> aVertex5(+0.000E+000, +7.500E-001, +0.000E+000);
    CGeoCoordinate<decimal> aVertex6(+0.000E+000, +7.500E-001, +2.500E-001);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-12);

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

    aTriangle2.intersects(aTriangle1, anIntersectionType, 1E-12);

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

}

TEST_F(CTestGeoTriangle, intersect9) {

    CGeoCoordinate<decimal> aVertex1(+0.000E+000, +5.000E-001, +2.500E-001);
    CGeoCoordinate<decimal> aVertex2(-6.100E-002, +6.667E-001, +2.233E-002);
    CGeoCoordinate<decimal> aVertex3(+2.500E-001, +7.500E-001, +0.000E+000);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(+0.000E+000, +5.000E-001, +0.000E+000);
    CGeoCoordinate<decimal> aVertex5(+0.000E+000, +7.500E-001, +0.000E+000);
    CGeoCoordinate<decimal> aVertex6(+0.000E+000, +7.500E-001, +2.500E-001);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-12);

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

    aTriangle2.intersects(aTriangle1, anIntersectionType, 1E-12);

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

}

TEST_F(CTestGeoTriangle, intersect10) {

    CGeoCoordinate<decimal> aVertex1(+7.375, -12.800, -8.499);
    CGeoCoordinate<decimal> aVertex2(+8.225, -12.670, -7.275);
    CGeoCoordinate<decimal> aVertex3(+8.650, -13.730, -8.500);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(+8.649, -12.96, -8.500);
    CGeoCoordinate<decimal> aVertex5(+8.650, -13.73, -8.500);
    CGeoCoordinate<decimal> aVertex6(+7.812, -13.68, -8.500);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-2);

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-12);

    EXPECT_TRUE(anIntersectionType == IT_VERTEX);

}

TEST_F(CTestGeoTriangle, intersect11) {

    CGeoCoordinate<decimal> aVertex1(+9.503E+000, -1.374E+001, -8.500E+000);
    CGeoCoordinate<decimal> aVertex2(+8.665E+000, -1.373E+001, -8.492E+000);
    CGeoCoordinate<decimal> aVertex3(+8.930E+000, -1.349E+001, -9.395E+000);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CGeoCoordinate<decimal> aVertex4(+9.503E+000, -1.374E+001, -8.500E+000);
    CGeoCoordinate<decimal> aVertex5(+8.665E+000, -1.373E+001, -8.500E+000);
    CGeoCoordinate<decimal> aVertex6(+8.649E+000, -1.296E+001, -8.500E+000);

    CGeoTriangle<decimal> aTriangle2;

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    CGeoIntersectionType anIntersectionType;

    aTriangle1.intersects(aTriangle2, anIntersectionType, 1E-2);

    EXPECT_TRUE(anIntersectionType == IT_EDGE);

}

TEST_F(CTestGeoTriangle, distance1) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 1.0, 0.0);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    decimal dist = 1.0;

    CGeoCoordinate<decimal> aPoint(0.1, 0.1, 0.0);
    CGeoCoordinate<decimal> aNewPoint;

    aTriangle1.distance(aPoint, aNewPoint, dist, 1E-6);

    EXPECT_LT(dist, 1E-6);

}

TEST_F(CTestGeoTriangle, distance2) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 1.0, 0.0);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    decimal dist = 0.0;

    CGeoCoordinate<decimal> aPoint(0.1, 0.1, 0.1);
    CGeoCoordinate<decimal> aNewPoint;

    aTriangle1.distance(aPoint, aNewPoint, dist, 1E-6);

    EXPECT_NEAR(dist, 0.1, 1E-5);

}

TEST_F(CTestGeoTriangle, distance3) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 1.0, 1.0);

    CGeoTriangle<decimal> aTriangle1;

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    decimal dist = 0.0;

    CGeoCoordinate<decimal> aPoint(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aNewPoint;

    aTriangle1.distance(aPoint, aNewPoint, dist, 1E-6);

    EXPECT_GT(dist, 1.0);

}


