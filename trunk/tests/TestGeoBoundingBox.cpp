// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#include "gtest/gtest.h"

#include "TypeDef.hpp"

#include "GeoBoundingBox.hpp"

using namespace ENigMA::geometry;

class CTestGeoBoundingBox : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoBoundingBox, contains) {

    CGeoBoundingBox<decimal> aBoundingBox;

    CGeoCoordinate<decimal> aCoordinate1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate2(1.0, 1.0, 1.0);
    
    aBoundingBox.addCoordinate(aCoordinate1);
    aBoundingBox.addCoordinate(aCoordinate2);

    CGeoCoordinate<decimal> aCoordinate3(-0.1, 0.0, 0.0);

    EXPECT_FALSE(aBoundingBox.contains(aCoordinate3));

}

TEST_F(CTestGeoBoundingBox, intersects) {

    CGeoBoundingBox<decimal> aBoundingBox1;

    CGeoCoordinate<decimal> aCoordinate1(+1.0, +6.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate2(-1.0, +2.0, 0.0);

    aBoundingBox1.addCoordinate(aCoordinate1);
    aBoundingBox1.addCoordinate(aCoordinate2);

    CGeoBoundingBox<decimal> aBoundingBox2;

    CGeoCoordinate<decimal> aCoordinate3(+1.0, +0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate4(-1.0, +4.0, 0.0);

    aBoundingBox2.addCoordinate(aCoordinate3);
    aBoundingBox2.addCoordinate(aCoordinate4);

    EXPECT_TRUE(aBoundingBox1.intersects(aBoundingBox2, 0.0));

}
