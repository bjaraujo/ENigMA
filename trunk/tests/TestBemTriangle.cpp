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

#include "BemTriangle.hpp"

using namespace ENigMA::bem;

class CTestBemTriangle : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestBemTriangle, integPoints) {

    CBemTriangle<Decimal> aTriangle;

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(0.0, 1.0, 1.0);

    aTriangle.addVertex(aVertex1);
    aTriangle.addVertex(aVertex2);
    aTriangle.addVertex(aVertex3);

    std::vector<CGeoCoordinate<Decimal> > sPoints;
    std::vector<Decimal> sWeights;

    sPoints.clear();
    sWeights.clear();

    aTriangle.getIntegrationPoints(sPoints, sWeights);

    EXPECT_EQ(20, sPoints.size());
    EXPECT_EQ(20, sWeights.size());

}

TEST_F(CTestBemTriangle, laplacianCoeff) {


    CBemTriangle<Decimal> aTriangle1;

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(0.0, 1.0, 0.0);

    aTriangle1.addVertex(aVertex1);
    aTriangle1.addVertex(aVertex2);
    aTriangle1.addVertex(aVertex3);

    CBemTriangle<Decimal> aTriangle2;

    CGeoCoordinate<Decimal> aVertex4(2.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(3.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex6(0.0, 1.0, 0.0);

    aTriangle2.addVertex(aVertex4);
    aTriangle2.addVertex(aVertex5);
    aTriangle2.addVertex(aVertex6);

    Decimal h, g;

    aTriangle1.laplacianCoeff(1, 2, aTriangle2, h, g);

    EXPECT_EQ(0.0, h);
    EXPECT_NE(0.0, g);

}

