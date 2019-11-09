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

#include "GeoVector.hpp"

using namespace ENigMA::geometry;

class CTestGeoVector : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoVector, set) {

    CGeoVector<decimal> aVector(-1.1, +2.2, -3.3);

    EXPECT_NEAR(-1.1, aVector.x(), 1E-6);
    EXPECT_NEAR(+2.2, aVector.y(), 1E-6);
    EXPECT_NEAR(-3.3, aVector.z(), 1E-6);

}

TEST_F(CTestGeoVector, operators) {

    CGeoVector<decimal> aVector1(-1.1, +2.2, -3.3);
    CGeoVector<decimal> aVector2(-1.1, +2.2, -3.3);

    EXPECT_TRUE(aVector1 == aVector2);

}

TEST_F(CTestGeoVector, angle) {

    CGeoVector<decimal> aVector1(+1.0, +0.0, +0.0);
    CGeoVector<decimal> aVector2(+0.0, +1.0, +0.0);

    const decimal pi = std::acos(-1.0);

    EXPECT_EQ(90, floor(aVector1.angle(aVector2) * 180.0 / pi + 0.5));

}

TEST_F(CTestGeoVector, cross) {

    CGeoVector<decimal> aVector1(+1.0, +0.0, +0.0);
    CGeoVector<decimal> aVector2(+0.0, +1.0, +0.0);
    CGeoVector<decimal> aVector3;

    aVector3 = aVector1.cross(aVector2);

    EXPECT_EQ(+0.0, aVector3.x());
    EXPECT_EQ(+0.0, aVector3.y());
    EXPECT_EQ(+1.0, aVector3.z());

}

TEST_F(CTestGeoVector, dot) {

    CGeoVector<decimal> aVector1(+1.0, +0.0, +0.0);
    CGeoVector<decimal> aVector2(-1.0, +0.0, +0.0);
    CGeoVector<decimal> aVector3;

    EXPECT_EQ(-1.0,  aVector1.dot(aVector2));

}

TEST_F(CTestGeoVector, rotate) {

    CGeoVector<decimal> aVector1(+1.0, +0.0, +0.0);

    const decimal pi = std::acos(-1.0);

    aVector1.rotate(pi/2.0);

    EXPECT_NEAR(0.0, aVector1.x(), 1E-6);

}
