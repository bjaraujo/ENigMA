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
#include "MshElement.hpp"

using namespace ENigMA::mesh;

class CTestMshElement : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestMshElement, set) {

    CMshElement<decimal> anElement;

    anElement.addNodeId(2);
    anElement.addNodeId(4);
    anElement.addNodeId(1);

    //std::cout << anElement << std::endl;

    EXPECT_EQ(3, anElement.nbNodeIds());

    EXPECT_EQ(2, anElement.nodeId(0));
    EXPECT_EQ(4, anElement.nodeId(1));
    EXPECT_EQ(1, anElement.nodeId(2));

}

