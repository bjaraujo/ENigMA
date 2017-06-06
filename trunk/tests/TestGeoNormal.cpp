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

#include "GeoNormal.hpp"

using namespace ENigMA::geometry;

class CTestGeoNormal : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoNormal, normalize) {

    CGeoNormal<decimal> aNormal(2.0,0.0,0.0);

    EXPECT_EQ(1.0, aNormal.x());
    EXPECT_EQ(0.0, aNormal.y());
    EXPECT_EQ(0.0, aNormal.z());

}

