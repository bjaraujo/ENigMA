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

#include "FemTetrahedron.hpp"

using namespace ENigMA::fem;

class CTestFemTetrahedron : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFemTetrahedron, update) {

    CFemTetrahedron<decimal, 4, 1, 1> aTetrahedron;

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.0, 1.0);

    aTetrahedron.addVertex(aVertex1);
    aTetrahedron.addVertex(aVertex2);
    aTetrahedron.addVertex(aVertex3);
    aTetrahedron.addVertex(aVertex4);

    aTetrahedron.calculateVolume();

    EXPECT_NEAR(1.0 / 6.0, aTetrahedron.volume(), 1E-6);

    aTetrahedron.transient() = true;

    aTetrahedron.dt() = 0.1;

    aTetrahedron.update();

    //std::cout << aTetrahedron.laplacian << std::endl;

    // TODO:
    EXPECT_NEAR(-1.0/3.0, aTetrahedron.laplacian(0, 0), 1E-6);

    //std::cout << aTetrahedron.ddt << std::endl;

}


