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

#include "GeoLine.hpp"

using namespace ENigMA::geometry;

class CTestGeoLine : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoLine, getLength) {

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(2.5, 0.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    aLine.calculateLength();

    EXPECT_EQ(2.5, aLine.length());

}

TEST_F(CTestGeoLine, clip1) {

    CGeoNormal<Decimal> aNormal(1.0, 0.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, 0.2);

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(0.0, aNewLine.startPoint().x(), 1E-6);
    EXPECT_NEAR(0.2, aNewLine.endPoint().x(), 1E-6);

    EXPECT_NEAR(0.2, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, clip2) {

    CGeoNormal<Decimal> aNormal(1.0, 0.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, 0.2);

    CGeoCoordinate<Decimal> aPoint1(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 1.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_EQ(0.0, aNewLine.length());

}

TEST_F(CTestGeoLine, clip3) {

    CGeoNormal<Decimal> aNormal(1.0, 0.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, 0.2);

    CGeoCoordinate<Decimal> aPoint1(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(0.0, 1.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(0.2, aNewLine.startPoint().x(), 1E-6);
    EXPECT_NEAR(0.0, aNewLine.endPoint().x(), 1E-6);

    EXPECT_NEAR(0.2, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, clip4) {

    CGeoNormal<Decimal> aNormal(1.0, 0.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, 0.2);

    CGeoCoordinate<Decimal> aPoint1(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(0.0, 0.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(1.0, aNewLine.startPoint().y(), 1E-6);
    EXPECT_NEAR(0.0, aNewLine.endPoint().y(), 1E-6);

    EXPECT_EQ(1.0, aNewLine.length());

}

TEST_F(CTestGeoLine, clip5) {

    CGeoNormal<Decimal> aNormal(1.0, 1.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, sqrt(2.0) / 2.0);

    CGeoCoordinate<Decimal> aPoint1(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 1.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(0.0, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, clip6) {

    CGeoNormal<Decimal> aNormal(1.0, 1.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, sqrt(2.0) / 2.0);

    CGeoCoordinate<Decimal> aPoint1(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(0.0, 1.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(0.0, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, clip7) {

    CGeoNormal<Decimal> aNormal(1.0, 1.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, sqrt(2.0) / 2.0);

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(0.0, aNewLine.startPoint().x(), 1E-6);
    EXPECT_NEAR(1.0, aNewLine.endPoint().x(), 1E-6);

    EXPECT_EQ(1.0, aNewLine.length());

}

TEST_F(CTestGeoLine, clip8) {

    CGeoNormal<Decimal> aNormal(1.0, 0.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, 1.5);

    CGeoCoordinate<Decimal> aPoint1(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(2.0, 0.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(1.0, aNewLine.startPoint().x(), 1E-6);
    EXPECT_NEAR(1.5, aNewLine.endPoint().x(), 1E-6);

    EXPECT_NEAR(0.5, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, clip9) {

    CGeoNormal<Decimal> aNormal(1.0, 1.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, 3.0 * sqrt(2.0) / 4.0);

    CGeoCoordinate<Decimal> aPoint1(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 1.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(0.0, aNewLine.startPoint().y(), 1E-6);
    EXPECT_NEAR(0.5, aNewLine.endPoint().y(), 1E-6);

    EXPECT_NEAR(0.5, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, clip10) {

    CGeoNormal<Decimal> aNormal(1.0, 1.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, 2.0);

    CGeoCoordinate<Decimal> aPoint1(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(0.0, 1.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(1.0, aNewLine.startPoint().x(), 1E-6);
    EXPECT_NEAR(0.0, aNewLine.endPoint().x(), 1E-6);

    EXPECT_NEAR(1.0, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, clip11) {

    CGeoNormal<Decimal> aNormal(0.0, 0.0, 1.0);
    CGeoPlane<Decimal> aPlane(aNormal, 0.5);

    CGeoCoordinate<Decimal> aPoint1(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 1.0, 1.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(0.0, aNewLine.startPoint().z(), 1E-6);
    EXPECT_NEAR(0.5, aNewLine.endPoint().z(), 1E-6);

    EXPECT_NEAR(0.5, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, clip12) {

    CGeoNormal<Decimal> aNormal(1.0, 0.0, 0.0);
    CGeoPlane<Decimal> aPlane(aNormal, 1.3);

    CGeoCoordinate<Decimal> aPoint1(2.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(0.0, 1.0, 0.0);

    CGeoLine<Decimal> aLine;

    aLine.setStartPoint(aPoint1);
    aLine.setEndPoint(aPoint2);

    CGeoLine<Decimal> aNewLine = aLine.clip(aPlane, 1E-9);

    aNewLine.calculateLength();

    EXPECT_NEAR(1.3, aNewLine.startPoint().x(), 1E-6);
    EXPECT_NEAR(0.0, aNewLine.endPoint().x(), 1E-6);

    EXPECT_NEAR(1.3, aNewLine.length(), 1E-6);

}

TEST_F(CTestGeoLine, intersection1) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(+1.0, +6.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(-1.0, +2.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(+1.0, +0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(-1.0, +4.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_NEAR(-0.5, aPoint.x(), 1E-6);
    EXPECT_NEAR(+3.0, aPoint.y(), 1E-6);

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

}

TEST_F(CTestGeoLine, intersection2) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(+1.0, +6.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(-1.0, +2.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(+1.0, +6.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(-1.0, +2.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_COINCIDENT);

}

TEST_F(CTestGeoLine, intersection3) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(+1.0, +6.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(-1.0, +2.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(+1.0, +6.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(-2.0, +2.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_VERTEX);

}

TEST_F(CTestGeoLine, intersection4) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(0.5, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(0.5, 1.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_FALSE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));
    
    EXPECT_TRUE(anIntersectionType == IT_NONE);

}

TEST_F(CTestGeoLine, intersection5) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(0.5, 0.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

}

TEST_F(CTestGeoLine, intersection6) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(0.1, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(0.5, 0.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

}

TEST_F(CTestGeoLine, intersection7) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(+0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(+1.0, 0.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(+0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(-1.0, 0.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_VERTEX);

}

TEST_F(CTestGeoLine, intersection8) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(+0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(+1.0, 0.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(+1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(+2.0, 0.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_VERTEX);

}

TEST_F(CTestGeoLine, intersection9) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(+0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(+1.0, 0.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(+0.5, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(+2.0, 0.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

}

TEST_F(CTestGeoLine, intersection10) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(+0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(+1.0, 0.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(-1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(+2.0, 0.0, 0.0);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_TRUE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_INTERNAL);

}

TEST_F(CTestGeoLine, intersection11) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(0, 0.125, 0);
    CGeoCoordinate<Decimal> aPoint2(0, 0, 0.125);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(0, 0, 0);
    CGeoCoordinate<Decimal> aPoint4(0, 0.125, -0.125);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint;

    CGeoIntersectionType anIntersectionType;

    EXPECT_FALSE(aLine1.intersects(aLine2, aPoint, anIntersectionType, 0.0));

    EXPECT_TRUE(anIntersectionType == IT_NONE);

}

TEST_F(CTestGeoLine, distance1) {

    CGeoLine<Decimal> aLine1;

    CGeoCoordinate<Decimal> aPoint1(-1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(+1.0, 0.0, 0.0);

    aLine1.setStartPoint(aPoint1);
    aLine1.setEndPoint(aPoint2);

    CGeoLine<Decimal> aLine2;

    CGeoCoordinate<Decimal> aPoint3(0.0, -1.0, 0.1);
    CGeoCoordinate<Decimal> aPoint4(0.0, +1.0, 0.1);

    aLine2.setStartPoint(aPoint3);
    aLine2.setEndPoint(aPoint4);

    CGeoCoordinate<Decimal> aPoint5;
    CGeoCoordinate<Decimal> aPoint6;

    Decimal dist = 0.0;

    EXPECT_TRUE(aLine1.distance(aLine2, aPoint5, aPoint6, dist, 1E-5));

    EXPECT_NEAR(0.1, dist, 1E-6);

}

TEST_F(CTestGeoLine, distance2) {

    CGeoCoordinate<Decimal> aPoint1(10.0, 150.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(400.0, 400.0, 0.0);

    CGeoLine<Decimal> aLine1(aPoint1, aPoint2);

    CGeoCoordinate<Decimal> aPoint3(10.0, 10.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(400.0, 255.0, 0.0);

    CGeoLine<Decimal> aLine2(aPoint3, aPoint4);

    CGeoCoordinate<Decimal> aPoint5;
    CGeoCoordinate<Decimal> aPoint6;

    Decimal dist = 0.0;

    EXPECT_TRUE(aLine1.distance(aLine2, aPoint5, aPoint6, dist, 1E-6));

    EXPECT_NEAR(118.54867793223777, dist, 1E-6);

}

TEST_F(CTestGeoLine, distance3) {

    CGeoCoordinate<Decimal> aPoint1(10.0, 160.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(400.0, 400.0, 0.0);

    CGeoLine<Decimal> aLine1(aPoint1, aPoint2);

    CGeoCoordinate<Decimal> aPoint3(10.0, 10.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(400.0, 272.5, 0.0);

    CGeoLine<Decimal> aLine2(aPoint3, aPoint4);

    CGeoCoordinate<Decimal> aPoint5;
    CGeoCoordinate<Decimal> aPoint6;

    Decimal dist = 0.0;

    EXPECT_TRUE(aLine1.distance(aLine2, aPoint5, aPoint6, dist, 1E-6));

    EXPECT_LT(dist, 127);

}
