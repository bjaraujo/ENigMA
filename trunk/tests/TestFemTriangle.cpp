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

#include "FemTriangle.hpp"

using namespace ENigMA::fem;

class CTestFemTriangle : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFemTriangle, update) {

    CFemTriangle<Decimal, 3, 1, 1> aTriangle;

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(0.0, 1.0, 1.0);

    aTriangle.addVertex(aVertex1);
    aTriangle.addVertex(aVertex2);
    aTriangle.addVertex(aVertex3);

    aTriangle.calculateArea();

    aTriangle.setTransient(true);

    aTriangle.setDt(0.1);

    aTriangle.update();

    //std::cout << aTriangle.laplacian << std::endl;

    EXPECT_NEAR(-0.5, aTriangle.laplacian(0, 0), 1E-6);
    EXPECT_NEAR(+0.5, aTriangle.laplacian(0, 1), 1E-6);
    EXPECT_NEAR(-0.0, aTriangle.laplacian(0, 2), 1E-6);

    EXPECT_NEAR(+0.5, aTriangle.laplacian(1, 0), 1E-6);
    EXPECT_NEAR(-1.0, aTriangle.laplacian(1, 1), 1E-6);
    EXPECT_NEAR(+0.5, aTriangle.laplacian(1, 2), 1E-6);

    EXPECT_NEAR(-0.0, aTriangle.laplacian(2, 0), 1E-6);
    EXPECT_NEAR(+0.5, aTriangle.laplacian(2, 1), 1E-6);
    EXPECT_NEAR(-0.5, aTriangle.laplacian(2, 2), 1E-6);

    //std::cout << aTriangle.ddt << std::endl;

}

