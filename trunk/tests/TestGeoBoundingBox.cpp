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

    CGeoBoundingBox<Decimal> aBoundingBox;

    CGeoCoordinate<Decimal> aCoordinate1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aCoordinate2(1.0, 1.0, 1.0);
    
    aBoundingBox.addCoordinate(aCoordinate1);
    aBoundingBox.addCoordinate(aCoordinate2);

    CGeoCoordinate<Decimal> aCoordinate3(-0.1, 0.0, 0.0);

    EXPECT_FALSE(aBoundingBox.contains(aCoordinate3));

}

TEST_F(CTestGeoBoundingBox, intersects) {

    CGeoBoundingBox<Decimal> aBoundingBox1;

    CGeoCoordinate<Decimal> aCoordinate1(+1.0, +6.0, 0.0);
    CGeoCoordinate<Decimal> aCoordinate2(-1.0, +2.0, 0.0);

    aBoundingBox1.addCoordinate(aCoordinate1);
    aBoundingBox1.addCoordinate(aCoordinate2);

    CGeoBoundingBox<Decimal> aBoundingBox2;

    CGeoCoordinate<Decimal> aCoordinate3(+1.0, +0.0, 0.0);
    CGeoCoordinate<Decimal> aCoordinate4(-1.0, +4.0, 0.0);

    aBoundingBox2.addCoordinate(aCoordinate3);
    aBoundingBox2.addCoordinate(aCoordinate4);

    EXPECT_TRUE(aBoundingBox1.intersects(aBoundingBox2, 0.0));

}
