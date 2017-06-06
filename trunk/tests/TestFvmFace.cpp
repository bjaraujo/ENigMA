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

#include "FvmFace.hpp"

using namespace ENigMA::fvm;

class CTestFvmFace : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFvmFace, calcArea) {

    CFvmFace<decimal> aFace;

    CFvmNode<decimal> aNode1(0, 0, 0);
    CFvmNode<decimal> aNode2(1, 0, 0);
    CFvmNode<decimal> aNode3(1, 1, 0);
    CFvmNode<decimal> aNode4(0, 1, 0);

    aFace.addNode(aNode1);
    aFace.addNode(aNode2);
    aFace.addNode(aNode3);
    aFace.addNode(aNode4);

    aFace.close();

    aFace.calculateArea();

    EXPECT_NEAR(1.0, aFace.area(), 1E-12);

}


