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

#include "MshNode.hpp"

using namespace ENigMA::mesh;

class CTestMshNode : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestMshNode, set) {

    CMshNode<Decimal> aNode(-1.1,+2.2,+3.3);

    //std::cout << aNode << std::endl;

    EXPECT_NEAR(-1.1, aNode.x(), 1E-6);
    EXPECT_NEAR(+2.2, aNode.y(), 1E-6);
    EXPECT_NEAR(+3.3, aNode.z(), 1E-6);

}

