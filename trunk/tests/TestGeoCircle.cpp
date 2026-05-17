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

#include "GeoCircle.hpp"

using namespace ENigMA::geometry;

const Decimal pi = static_cast<Decimal>(std::acos(-1.0));

class CTestGeoCircle : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoCircle, circle1) {

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint3(0.5, 0.5, 0.0);

    CGeoCircle<Decimal> aCircle(aPoint1, aPoint2, aPoint3);

    aCircle.calculateArea();

    EXPECT_NEAR(pi/4.0, aCircle.area(), 1E-3);

}

TEST_F(CTestGeoCircle, boundingBox1) {

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint3(0.5, 0.5, 0.0);

    CGeoCircle<Decimal> aCircle(aPoint1, aPoint2, aPoint3);

    aCircle.calculateBoundingBox();

    EXPECT_NEAR(+0.0, aCircle.boundingBox().min().x(), 1E-3);
    EXPECT_NEAR(-0.5, aCircle.boundingBox().min().y(), 1E-3);
    EXPECT_NEAR(+0.0, aCircle.boundingBox().min().z(), 1E-3);

    EXPECT_NEAR(+1.0, aCircle.boundingBox().max().x(), 1E-3);
    EXPECT_NEAR(+0.5, aCircle.boundingBox().max().y(), 1E-3);
    EXPECT_NEAR(+0.0, aCircle.boundingBox().max().z(), 1E-3);
}

TEST_F(CTestGeoCircle, area1) {

    CGeoCoordinate<Decimal> aCenter(0.0, 0.0, 0.0);

    CGeoCircle<Decimal> aCircle(aCenter, 1.0);

    aCircle.calculateArea();

    EXPECT_NEAR(pi, aCircle.area(), 1E-3);

}

TEST_F(CTestGeoCircle, area2) {

    CGeoCoordinate<Decimal> aPoint1(+1.0, +0.0, +0.0);
    CGeoCoordinate<Decimal> aPoint2(+0.0, +1.0, +0.0);
    CGeoCoordinate<Decimal> aPoint3(-1.0, +0.0, +0.0);

    CGeoCircle<Decimal> aCircle(aPoint1, aPoint2, aPoint3);

    aCircle.calculateArea();

    EXPECT_NEAR(pi, aCircle.area(), 1E-3);

}

TEST_F(CTestGeoCircle, contains1) {

    CGeoCoordinate<Decimal> aPoint1(+1.0, +0.0, +0.0);
    CGeoCoordinate<Decimal> aPoint2(+0.0, +1.0, +0.0);
    CGeoCoordinate<Decimal> aPoint3(-1.0, +0.0, +0.0);

    CGeoCircle<Decimal> aCircle(aPoint1, aPoint2, aPoint3);

    CGeoCoordinate<Decimal> aPoint5(0.25, 0.25, 0.0);
    CGeoCoordinate<Decimal> aPoint6(1.5, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint7(1.5, 1.5, 1.5);

    EXPECT_TRUE(aCircle.contains(aPoint5));
    EXPECT_FALSE(aCircle.contains(aPoint6));
    EXPECT_FALSE(aCircle.contains(aPoint7));

}
