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

#include "FemHexahedron.hpp"

using namespace ENigMA::fem;

class CTestFemHexahedron : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFemHexahedron, update) {

    CFemHexahedron<decimal, 8, 1, 1> aHexahedron;

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex6(1.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex7(1.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aVertex8(0.0, 1.0, 1.0);

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    aHexahedron.calculateVolume();

    aHexahedron.setTransient(true);

    aHexahedron.setDt(0.1);

    aHexahedron.update();

    //std::cout << aHexahedron.laplacian << std::endl;
    //std::cout << aHexahedron.ddt << std::endl;

    EXPECT_NEAR(-0.352825164509670, aHexahedron.laplacian(0,0), 1E-6);
    EXPECT_NEAR(+0.160375074774608, aHexahedron.laplacian(1,0), 1E-6);
    EXPECT_NEAR(-0.128300059817612, aHexahedron.laplacian(2,0), 1E-6);

    EXPECT_NEAR(+0.213833433042801, aHexahedron.ddt(0,0), 1E-6);
    EXPECT_NEAR(+0.106916716518519, aHexahedron.ddt(1,0), 1E-6);
    EXPECT_NEAR(+0.106916716518519, aHexahedron.ddt(2,0), 1E-6);

}


