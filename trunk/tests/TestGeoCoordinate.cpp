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

#include "GeoCoordinate.hpp"

using namespace ENigMA::geometry;

class CTestGeoCoordinate : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoCoordinate, set1) {

    CGeoCoordinate<Decimal> aCoordinate;

    EXPECT_NEAR(+0.0, aCoordinate.x(), 1E-20);
    EXPECT_NEAR(+0.0, aCoordinate.y(), 1E-20);
    EXPECT_NEAR(+0.0, aCoordinate.z(), 1E-20);

}

TEST_F(CTestGeoCoordinate, set2) {

    CGeoCoordinate<Decimal> aCoordinate(-1.1, +2.2, -3.3);

    EXPECT_NEAR(-1.1, aCoordinate.x(), 1E-6);
    EXPECT_NEAR(+2.2, aCoordinate.y(), 1E-6);
    EXPECT_NEAR(-3.3, aCoordinate.z(), 1E-6);

}

TEST_F(CTestGeoCoordinate, transform) {

    CGeoCoordinate<Decimal> aCoordinate(-1.1, +2.2, -3.3);

    CGeoCoordinateSystem<Decimal> aCoordinateSystem;
    aCoordinateSystem << 1, 0, 0, 0, 1, 0, 0, 0, 1;

    aCoordinate.transform(aCoordinateSystem);

    EXPECT_NEAR(-1.1, aCoordinate.x(), 1E-6);
    EXPECT_NEAR(+2.2, aCoordinate.y(), 1E-6);
    EXPECT_NEAR(-3.3, aCoordinate.z(), 1E-6);

}
