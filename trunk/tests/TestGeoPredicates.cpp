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

#include "GeoPredicates.hpp"

using namespace ENigMA::geometry;

class CTestGeoPredicates : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoPredicates, inCircle) {

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint3(1.0, 1.0, 0.0);

    CGeoPredicates<decimal>::exactInit();

    CGeoCoordinate<decimal> aPoint(0.5, 0.5, 0.0);

    EXPECT_GT(CGeoPredicates<decimal>::inCircle(aPoint1, aPoint2, aPoint3, aPoint), 0.0);

}

TEST_F(CTestGeoPredicates, inSphere) {

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aPoint3(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint4(0.0, 0.0, 1.0);
    
    CGeoCoordinate<decimal> aPoint;
    
    CGeoPredicates<decimal>::exactInit();
    
    aPoint << 0.5, 0.5, 0.5;

    EXPECT_GT(CGeoPredicates<decimal>::inSphere(aPoint1, aPoint2, aPoint3, aPoint4, aPoint), 0.0);

}

