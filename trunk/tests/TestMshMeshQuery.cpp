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

#include "MshMeshQuery.hpp"

using namespace ENigMA::mesh;

class CTestMshMeshQuery : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestMshMeshQuery, query) {

    CMshMesh<decimal> aMesh;

    CMshNode<decimal> aNode0(0,0,0);
    aMesh.addNode(0, aNode0);

    CMshNode<decimal> aNode1(1,0,0);
    aMesh.addNode(1, aNode1);

    CMshNode<decimal> aNode2(1,1,0);
    aMesh.addNode(2, aNode2);

    CMshNode<decimal> aNode3(0,1,0);
    aMesh.addNode(3, aNode3);

    CMshNode<decimal> aNode4(1,-1,0);
    aMesh.addNode(4, aNode4);

    CMshNode<decimal> aNode5(-1,-1,0);
    aMesh.addNode(5, aNode5);

    CMshNode<decimal> aNode6(-2,-1,0);
    aMesh.addNode(6, aNode6);
    
    CMshElement<decimal> anElement0(ET_TRIANGLE);
    anElement0.addNodeId(0);
    anElement0.addNodeId(1);
    anElement0.addNodeId(2);
    aMesh.addElement(0, anElement0);

    CMshElement<decimal> anElement1(ET_TRIANGLE);
    anElement1.addNodeId(0);
    anElement1.addNodeId(2);
    anElement1.addNodeId(3);
    aMesh.addElement(1, anElement1);

    CMshElement<decimal> anElement2(ET_TRIANGLE);
    anElement2.addNodeId(0);
    anElement2.addNodeId(4);
    anElement2.addNodeId(1);
    aMesh.addElement(2, anElement2);

    CMshElement<decimal> anElement3(ET_TRIANGLE);
    anElement3.addNodeId(0);
    anElement3.addNodeId(5);
    anElement3.addNodeId(4);
    aMesh.addElement(3, anElement3);

    CMshElement<decimal> anElement4(ET_TRIANGLE);
    anElement4.addNodeId(5);
    anElement4.addNodeId(6);
    anElement4.addNodeId(4);
    aMesh.addElement(4, anElement4);
    
    CMshMeshQuery<decimal> aMeshQuery(aMesh);
    
    std::vector<Integer> sElementIds;
    
    aMeshQuery.elementsSharingNode(0, sElementIds);
    
    EXPECT_EQ(4, sElementIds.size());
    
    EXPECT_TRUE(std::find(sElementIds.begin(), sElementIds.end(), 0) != sElementIds.end()); 
    EXPECT_TRUE(std::find(sElementIds.begin(), sElementIds.end(), 1) != sElementIds.end()); 
    EXPECT_TRUE(std::find(sElementIds.begin(), sElementIds.end(), 2) != sElementIds.end()); 
    EXPECT_TRUE(std::find(sElementIds.begin(), sElementIds.end(), 3) != sElementIds.end()); 

    aMeshQuery.elementsSharingNodes(0, 2, sElementIds);

    EXPECT_EQ(2, sElementIds.size());

    EXPECT_TRUE(std::find(sElementIds.begin(), sElementIds.end(), 0) != sElementIds.end()); 
    EXPECT_TRUE(std::find(sElementIds.begin(), sElementIds.end(), 1) != sElementIds.end()); 

}

