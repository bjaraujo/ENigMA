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

#include "FemTriangularPrism.hpp"

using namespace ENigMA::fem;

class CTestFemTriangularPrism : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFemTriangularPrism, update) {

    CFemTriangularPrism<Decimal, 6, 1, 1> aTriangularPrism;

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, -0.5);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, -0.5);
    CGeoCoordinate<Decimal> aVertex3(0.0, 1.0, -0.5);
    CGeoCoordinate<Decimal> aVertex4(0.0, 0.0, +0.5);
    CGeoCoordinate<Decimal> aVertex5(1.0, 0.0, +0.5);
    CGeoCoordinate<Decimal> aVertex6(0.0, 1.0, +0.5);

    aTriangularPrism.addVertex(aVertex1);
    aTriangularPrism.addVertex(aVertex2);
    aTriangularPrism.addVertex(aVertex3);
    aTriangularPrism.addVertex(aVertex4);
    aTriangularPrism.addVertex(aVertex5);
    aTriangularPrism.addVertex(aVertex6);

    aTriangularPrism.calculateVolume();

    aTriangularPrism.setTransient(true);

    aTriangularPrism.setDt(0.1);

    aTriangularPrism.update();

    //std::cout << aTriangularPrism.laplacian << std::endl;
    //std::cout << aTriangularPrism.ddt << std::endl;

    EXPECT_NEAR(-0.368055555555556, aTriangularPrism.laplacian(0,0), 1E-6);
    EXPECT_NEAR(+0.0798611111111111, aTriangularPrism.laplacian(1,0), 1E-6);
    EXPECT_NEAR(+0.0798611111111111, aTriangularPrism.laplacian(2,0), 1E-6);

    EXPECT_NEAR(+0.393518518518518, aTriangularPrism.ddt(0,0), 1E-6);
    EXPECT_NEAR(+0.150462962962963, aTriangularPrism.ddt(1,0), 1E-6);
    EXPECT_NEAR(+0.150462962962963, aTriangularPrism.ddt(2,0), 1E-6);

}


