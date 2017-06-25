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

#include "GeoSphere.hpp"

using namespace ENigMA::geometry;

class CTestGeoSphere : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoSphere, volume1) {

    CGeoCoordinate<decimal> aCenter(0.0, 0.0, 0.0);

    CGeoSphere<decimal> aSphere(aCenter, 1.0);

    aSphere.calculateVolume();

    EXPECT_NEAR(4.18879020479, aSphere.volume(), 1E-3);

}

TEST_F(CTestGeoSphere, volume2) {

    CGeoCoordinate<decimal> aPoint1(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aPoint2(+0.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aPoint3(-1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aPoint4(+0.0, +0.0, +1.0);

    CGeoSphere<decimal> aSphere(aPoint1, aPoint2, aPoint3, aPoint4);

    aSphere.calculateVolume();

    EXPECT_NEAR(4.18879020479, aSphere.volume(), 1E-3);

}

TEST_F(CTestGeoSphere, contains1) {

    CGeoCoordinate<decimal> aPoint1(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aPoint2(+0.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aPoint3(-1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aPoint4(+0.0, +0.0, +1.0);

    CGeoSphere<decimal> aSphere(aPoint1, aPoint2, aPoint3, aPoint4);

    CGeoCoordinate<decimal> aPoint5(0, 0, 0);
    CGeoCoordinate<decimal> aPoint6(0, 0, 1.5);
    CGeoCoordinate<decimal> aPoint7(1.5, 1.5, 1.5);

    EXPECT_TRUE(aSphere.contains(aPoint5));
    EXPECT_FALSE(aSphere.contains(aPoint6));
    EXPECT_FALSE(aSphere.contains(aPoint7));

}
