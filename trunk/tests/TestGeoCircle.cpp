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

class CTestGeoCircle : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoCircle, area1) {

    CGeoCoordinate<decimal> aCenter(0.0, 0.0, 0.0);

    CGeoCircle<decimal> aCircle(aCenter, 1.0);

    aCircle.calculateArea();

    EXPECT_NEAR(3.14159265359, aCircle.area(), 1E-3);

}

TEST_F(CTestGeoCircle, area2) {

    CGeoCoordinate<decimal> aPoint1(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aPoint2(+0.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aPoint3(-1.0, +0.0, +0.0);

    CGeoCircle<decimal> aCircle(aPoint1, aPoint2, aPoint3);

    aCircle.calculateArea();

    EXPECT_NEAR(3.14159265359, aCircle.area(), 1E-3);

}

TEST_F(CTestGeoCircle, contains1) {

    CGeoCoordinate<decimal> aPoint1(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aPoint2(+0.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aPoint3(-1.0, +0.0, +0.0);

    CGeoCircle<decimal> aCircle(aPoint1, aPoint2, aPoint3);

    CGeoCoordinate<decimal> aPoint5(0.25, 0.25, 0.0);
    CGeoCoordinate<decimal> aPoint6(1.5, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint7(1.5, 1.5, 1.5);

    EXPECT_TRUE(aCircle.contains(aPoint5));
    EXPECT_FALSE(aCircle.contains(aPoint6));
    EXPECT_FALSE(aCircle.contains(aPoint7));

}
