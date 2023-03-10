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

#include "GeoPlane.hpp"

using namespace ENigMA::geometry;

class CTestGeoPlane : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoPlane, construct) {

    CGeoNormal<Decimal> aNormal(2, 2, 2);

    CGeoPlane<Decimal> aPlane;

    aPlane.set(aNormal, 1.0);

    EXPECT_NEAR(1.0, aPlane.normal().norm(), 1E-6);
    EXPECT_NEAR(1.0, aPlane.d(), 1E-6);

}

TEST_F(CTestGeoPlane, set) {


    CGeoNormal<Decimal> aNormal(2, 2, 2);

    CGeoPlane<Decimal> aPlane(aNormal, 0.0);

    EXPECT_NEAR(1.0, aPlane.normal().norm(), 1E-6);
    EXPECT_NEAR(0.0, aPlane.d(), 1E-6);

}

TEST_F(CTestGeoPlane, distance) {

    CGeoNormal<Decimal> aNormal(0.0, 0.0, 1.0);

    CGeoPlane<Decimal> aPlane(aNormal, 0.1);

    CGeoCoordinate<Decimal> aPoint(0.0, 0.0, 0.0);

    CGeoCoordinate<Decimal> aNewPoint;

    Decimal dist = aPlane.distance(aPoint, aNewPoint);

    EXPECT_NEAR(0.1, dist, 1E-6);

    EXPECT_NEAR(0.0, aNewPoint.x(), 1E-6);
    EXPECT_NEAR(0.0, aNewPoint.y(), 1E-6);
    EXPECT_NEAR(0.1, aNewPoint.z(), 1E-6);

}
