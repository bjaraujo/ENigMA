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

    CGeoCoordinate<Decimal> aCenter(0.0, 0.0, 0.0);

    CGeoSphere<Decimal> aSphere(aCenter, 1.0);

    aSphere.calculateVolume();

    EXPECT_NEAR(4.18879020479, aSphere.volume(), 1E-3);

}

TEST_F(CTestGeoSphere, volume2) {

    CGeoCoordinate<Decimal> aPoint1(+1.0, +0.0, +0.0);
    CGeoCoordinate<Decimal> aPoint2(+0.0, +1.0, +0.0);
    CGeoCoordinate<Decimal> aPoint3(-1.0, +0.0, +0.0);
    CGeoCoordinate<Decimal> aPoint4(+0.0, +0.0, +1.0);

    CGeoSphere<Decimal> aSphere(aPoint1, aPoint2, aPoint3, aPoint4);

    aSphere.calculateVolume();

    EXPECT_NEAR(4.18879020479, aSphere.volume(), 1E-3);

}

TEST_F(CTestGeoSphere, contains1) {

    CGeoCoordinate<Decimal> aPoint1(+1.0, +0.0, +0.0);
    CGeoCoordinate<Decimal> aPoint2(+0.0, +1.0, +0.0);
    CGeoCoordinate<Decimal> aPoint3(-1.0, +0.0, +0.0);
    CGeoCoordinate<Decimal> aPoint4(+0.0, +0.0, +1.0);

    CGeoSphere<Decimal> aSphere(aPoint1, aPoint2, aPoint3, aPoint4);

    CGeoCoordinate<Decimal> aPoint5(0, 0, 0);
    CGeoCoordinate<Decimal> aPoint6(0, 0, 1.5);
    CGeoCoordinate<Decimal> aPoint7(1.5, 1.5, 1.5);

    EXPECT_TRUE(aSphere.contains(aPoint5));
    EXPECT_FALSE(aSphere.contains(aPoint6));
    EXPECT_FALSE(aSphere.contains(aPoint7));

}
