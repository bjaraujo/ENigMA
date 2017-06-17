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

#include "FemQuadrilateral.hpp"

using namespace ENigMA::fem;

class CTestFemQuadrilateral : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFemQuadrilateral, update) {

    CFemQuadrilateral<decimal, 4, 1, 1> aQuadrilateral;

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.0, 1.0);

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    aQuadrilateral.calculateArea();

    aQuadrilateral.transient() = true;

    aQuadrilateral.dt() = 0.1;

    aQuadrilateral.update();

    //std::cout << aQuadrilateral.laplacian << std::endl;

    EXPECT_NEAR(-2.0/3.0, aQuadrilateral.laplacian(0, 0), 1E-6);
    EXPECT_NEAR(+1.0/6.0, aQuadrilateral.laplacian(0, 1), 1E-6);
    EXPECT_NEAR(+1.0/3.0, aQuadrilateral.laplacian(0, 2), 1E-6);
    EXPECT_NEAR(+1.0/6.0, aQuadrilateral.laplacian(0, 3), 1E-6);

    EXPECT_NEAR(+1.0/6.0, aQuadrilateral.laplacian(1, 0), 1E-6);
    EXPECT_NEAR(-2.0/3.0, aQuadrilateral.laplacian(1, 1), 1E-6);
    EXPECT_NEAR(+1.0/6.0, aQuadrilateral.laplacian(1, 2), 1E-6);
    EXPECT_NEAR(+1.0/3.0, aQuadrilateral.laplacian(1, 3), 1E-6);

    EXPECT_NEAR(+1.0/3.0, aQuadrilateral.laplacian(2, 0), 1E-6);
    EXPECT_NEAR(+1.0/6.0, aQuadrilateral.laplacian(2, 1), 1E-6);
    EXPECT_NEAR(-2.0/3.0, aQuadrilateral.laplacian(2, 2), 1E-6);
    EXPECT_NEAR(+1.0/6.0, aQuadrilateral.laplacian(2, 3), 1E-6);

    EXPECT_NEAR(+1.0/6.0, aQuadrilateral.laplacian(3, 0), 1E-6);
    EXPECT_NEAR(+1.0/3.0, aQuadrilateral.laplacian(3, 1), 1E-6);
    EXPECT_NEAR(+1.0/6.0, aQuadrilateral.laplacian(3, 2), 1E-6);
    EXPECT_NEAR(-2.0/3.0, aQuadrilateral.laplacian(3, 3), 1E-6);

    //std::cout << aQuadrilateral.ddt << std::endl;

}


