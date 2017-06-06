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

#include "GeoTriangularPrism.hpp"

using namespace ENigMA::geometry;

class CTestGeoTriangularPrism : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoTriangularPrism, volume) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);

    CGeoCoordinate<decimal> aVertex4(0.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex5(1.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex6(1.0, 1.0, 1.0);

    CGeoTriangularPrism<decimal> aTriangularPrism;

    aTriangularPrism.addVertex(aVertex1);
    aTriangularPrism.addVertex(aVertex2);
    aTriangularPrism.addVertex(aVertex3);
    aTriangularPrism.addVertex(aVertex4);
    aTriangularPrism.addVertex(aVertex5);
    aTriangularPrism.addVertex(aVertex6);

    aTriangularPrism.calculateVolume();

    EXPECT_EQ(0.5, aTriangularPrism.volume());

}

