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

#include "FvmControlVolume.hpp"

using namespace ENigMA::fvm;

class CTestFvmControlVolume : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFvmControlVolume, calcVolume1) {

    CFvmControlVolume<decimal> aControlVolume;

    CFvmNode<decimal> aNode1(0.0, 0.0, 0.0);
    CFvmNode<decimal> aNode2(1.0, 0.0, 0.0);
    CFvmNode<decimal> aNode3(1.0, 1.0, 0.0);
    CFvmNode<decimal> aNode4(0.0, 1.0, 0.0);
    CFvmNode<decimal> aNode5(0.0, 0.0, 1.0);
    CFvmNode<decimal> aNode6(1.0, 0.0, 1.0);
    CFvmNode<decimal> aNode7(1.0, 1.0, 1.0);
    CFvmNode<decimal> aNode8(0.0, 1.0, 1.0);

    CFvmFace<decimal> aFace;

    // Face 1
    aFace.addNode(aNode1);
    aFace.addNode(aNode4);
    aFace.addNode(aNode3);
    aFace.addNode(aNode2);
    aFace.close();

    aControlVolume.addFace(0, aFace);
    aFace.reset();

    // Face 2
    aFace.addNode(aNode5);
    aFace.addNode(aNode6);
    aFace.addNode(aNode7);
    aFace.addNode(aNode8);
    aFace.close();

    aControlVolume.addFace(1, aFace);
    aFace.reset();

    // Face 3
    aFace.addNode(aNode5);
    aFace.addNode(aNode1);
    aFace.addNode(aNode2);
    aFace.addNode(aNode6);
    aFace.close();

    aControlVolume.addFace(2, aFace);
    aFace.reset();

    // Face 4
    aFace.addNode(aNode8);
    aFace.addNode(aNode7);
    aFace.addNode(aNode3);
    aFace.addNode(aNode4);
    aFace.close();

    aControlVolume.addFace(3, aFace);
    aFace.reset();

    // Face 5
    aFace.addNode(aNode6);
    aFace.addNode(aNode2);
    aFace.addNode(aNode3);
    aFace.addNode(aNode7);
    aFace.close();

    aControlVolume.addFace(4, aFace);
    aFace.reset();

    // Face 6
    aFace.addNode(aNode1);
    aFace.addNode(aNode5);
    aFace.addNode(aNode8);
    aFace.addNode(aNode4);
    aFace.close();

    aControlVolume.addFace(5, aFace);
    aFace.reset();

    aControlVolume.calculateVolume(true);

    EXPECT_NEAR(1.0, aControlVolume.volume(), 1E-12);

}

TEST_F(CTestFvmControlVolume, clip1) {

    CFvmControlVolume<decimal> aControlVolume;

    CFvmNode<decimal> aNode1(0.0, 0.0, 0.0);
    CFvmNode<decimal> aNode2(1.0, 0.0, 0.0);
    CFvmNode<decimal> aNode3(1.0, 1.0, 0.0);
    CFvmNode<decimal> aNode4(0.0, 1.0, 0.0);
    CFvmNode<decimal> aNode5(0.0, 0.0, 1.0);
    CFvmNode<decimal> aNode6(1.0, 0.0, 1.0);
    CFvmNode<decimal> aNode7(1.0, 1.0, 1.0);
    CFvmNode<decimal> aNode8(0.0, 1.0, 1.0);

    CFvmFace<decimal> aFace;

    // Face 1
    aFace.addNode(aNode1);
    aFace.addNode(aNode4);
    aFace.addNode(aNode3);
    aFace.addNode(aNode2);
    aFace.close();

    aControlVolume.addFace(0, aFace);
    aFace.reset();

    // Face 2
    aFace.addNode(aNode5);
    aFace.addNode(aNode6);
    aFace.addNode(aNode7);
    aFace.addNode(aNode8);
    aFace.close();

    aControlVolume.addFace(1, aFace);
    aFace.reset();

    // Face 3
    aFace.addNode(aNode5);
    aFace.addNode(aNode1);
    aFace.addNode(aNode2);
    aFace.addNode(aNode6);
    aFace.close();

    aControlVolume.addFace(2, aFace);
    aFace.reset();

    // Face 4
    aFace.addNode(aNode8);
    aFace.addNode(aNode7);
    aFace.addNode(aNode3);
    aFace.addNode(aNode4);
    aFace.close();

    aControlVolume.addFace(3, aFace);
    aFace.reset();

    // Face 5
    aFace.addNode(aNode6);
    aFace.addNode(aNode2);
    aFace.addNode(aNode3);
    aFace.addNode(aNode7);
    aFace.close();

    aControlVolume.addFace(4, aFace);
    aFace.reset();

    // Face 6
    aFace.addNode(aNode1);
    aFace.addNode(aNode5);
    aFace.addNode(aNode8);
    aFace.addNode(aNode4);
    aFace.close();

    aControlVolume.addFace(5, aFace);
    aFace.reset();

    aControlVolume.calculateVolume(true);

    EXPECT_NEAR(1.0, aControlVolume.volume(), 1E-12);

    CGeoNormal<decimal> aNormal(1.0, 0.3, 0.2);

    decimal volumeFractionAct;
    Integer nIterations;

    aControlVolume.setClippedFaceId(999);
    aControlVolume.clip(aNormal, 0.4, volumeFractionAct, nIterations, 50, 1E-6);

    aControlVolume.calculateVolume(true);

    EXPECT_NEAR(0.4, aControlVolume.volume(), 1E-3);

    aControlVolume.setClippedFaceId(999);
    aControlVolume.clip(aNormal, 0.8, volumeFractionAct, nIterations, 50, 1E-6);

    aControlVolume.calculateVolume(true);

    EXPECT_NEAR(0.8, aControlVolume.volume(), 1E-3);

}

TEST_F(CTestFvmControlVolume, clip2) {

    CFvmControlVolume<decimal> aControlVolume;

    CFvmNode<decimal> aNode1(0.0, 0.0, 0.0);
    CFvmNode<decimal> aNode2(1.0, 0.0, 0.0);
    CFvmNode<decimal> aNode3(1.0, 1.0, 0.0);
    CFvmNode<decimal> aNode4(0.0, 1.0, 0.0);
    CFvmNode<decimal> aNode5(0.0, 0.0, 1.0);
    CFvmNode<decimal> aNode6(1.0, 0.0, 1.0);
    CFvmNode<decimal> aNode7(1.0, 1.0, 1.0);
    CFvmNode<decimal> aNode8(0.0, 1.0, 1.0);

    CFvmFace<decimal> aFace;

    // Face 1
    aFace.addNode(aNode1);
    aFace.addNode(aNode4);
    aFace.addNode(aNode3);
    aFace.addNode(aNode2);
    aFace.close();

    aControlVolume.addFace(0, aFace);
    aFace.reset();

    // Face 2
    aFace.addNode(aNode5);
    aFace.addNode(aNode6);
    aFace.addNode(aNode7);
    aFace.addNode(aNode8);
    aFace.close();

    aControlVolume.addFace(1, aFace);
    aFace.reset();

    // Face 3
    aFace.addNode(aNode5);
    aFace.addNode(aNode1);
    aFace.addNode(aNode2);
    aFace.addNode(aNode6);
    aFace.close();

    aControlVolume.addFace(2, aFace);
    aFace.reset();

    // Face 4
    aFace.addNode(aNode8);
    aFace.addNode(aNode7);
    aFace.addNode(aNode3);
    aFace.addNode(aNode4);
    aFace.close();

    aControlVolume.addFace(3, aFace);
    aFace.reset();

    // Face 5
    aFace.addNode(aNode6);
    aFace.addNode(aNode2);
    aFace.addNode(aNode3);
    aFace.addNode(aNode7);
    aFace.close();
    aControlVolume.addFace(4, aFace);
    aFace.reset();

    // Face 6
    aFace.addNode(aNode1);
    aFace.addNode(aNode5);
    aFace.addNode(aNode8);
    aFace.addNode(aNode4);
    aFace.close();

    aControlVolume.addFace(5, aFace);
    aFace.reset();

    aControlVolume.calculateCentroid(true);

    EXPECT_NEAR(0.5, aControlVolume.centroid().x(), 1E-3);
    EXPECT_NEAR(0.5, aControlVolume.centroid().y(), 1E-3);
    EXPECT_NEAR(0.5, aControlVolume.centroid().z(), 1E-3);

    aControlVolume.calculateVolume(true);

    EXPECT_NEAR(1.0, aControlVolume.volume(), 1E-12);

    decimal sumArea1 = 0.0;
    decimal sumDistance1 = 0.0;

    for (Integer i = 0; i < aControlVolume.nbFaces(); ++i)
    {

        Integer aFaceId = aControlVolume.faceId(i);

        // Area
        aControlVolume.calculateFaceArea(aFaceId);

        sumArea1 += aControlVolume.faceArea(aFaceId);

        // Distance
        aControlVolume.calculateFaceCentroid(aFaceId);

        EXPECT_NEAR(0.5, aControlVolume.faceDist(aFaceId), 1E-3);

        sumDistance1 += aControlVolume.faceDist(aFaceId);

    }

    EXPECT_NEAR(6.0, sumArea1, 1E-3);
    EXPECT_NEAR(3.0, sumDistance1, 1E-3);

    CGeoNormal<decimal> aNormal(1.0, 0.0, 0.0);

    decimal volumeFractionAct;
    Integer nIterations;

    aControlVolume.setClippedFaceId(999);
    aControlVolume.clip(aNormal, 0.5, volumeFractionAct, nIterations, 50, 1E-6);

    aControlVolume.calculateCentroid(true);

    EXPECT_NEAR(0.25, aControlVolume.centroid().x(), 1E-3);
    EXPECT_NEAR(0.50, aControlVolume.centroid().y(), 1E-3);
    EXPECT_NEAR(0.50, aControlVolume.centroid().z(), 1E-3);

    aControlVolume.calculateVolume(true);

    EXPECT_NEAR(0.5, aControlVolume.volume(), 1E-3);

    decimal sumArea2 = 0.0;
    decimal sumNormal2 = 0.0;
    decimal sumDistance2 = 0.0;

    for (Integer i = 0; i < aControlVolume.nbFaces(); ++i)
    {

        Integer aFaceId = aControlVolume.faceId(i);

        // Area
        aControlVolume.calculateFaceArea(aFaceId);

        sumArea2 += aControlVolume.faceArea(aFaceId);

        // Normal
        sumNormal2 += aControlVolume.faceNormal(aFaceId).x() + 
            aControlVolume.faceNormal(aFaceId).y() + 
            aControlVolume.faceNormal(aFaceId).z();

        aControlVolume.calculateFaceCentroid(aFaceId);

        sumDistance2 += aControlVolume.faceDist(aFaceId);

    }

    EXPECT_NEAR(4.0, sumArea2, 1E-3);
    EXPECT_NEAR(2.5, sumDistance2, 1E-3);
    EXPECT_NEAR(0.0, sumNormal2, 1E-3);

}

TEST_F(CTestFvmControlVolume, clip3) {

    CFvmControlVolume<decimal> aControlVolume;

    CFvmNode<decimal> aNode1(0.0, 0.0, 0.0);
    CFvmNode<decimal> aNode2(1.0, 0.0, 0.0);
    CFvmNode<decimal> aNode3(1.0, 1.0, 0.0);
    CFvmNode<decimal> aNode4(0.0, 1.0, 0.0);
    CFvmNode<decimal> aNode5(0.0, 0.0, 1.0);
    CFvmNode<decimal> aNode6(1.0, 0.0, 1.0);
    CFvmNode<decimal> aNode7(1.0, 1.0, 1.0);
    CFvmNode<decimal> aNode8(0.0, 1.0, 1.0);

    CFvmFace<decimal> aFace;

    // Face 1
    aFace.addNode(aNode1);
    aFace.addNode(aNode4);
    aFace.addNode(aNode3);
    aFace.addNode(aNode2);
    aFace.close();

    aControlVolume.addFace(0, aFace);
    aFace.reset();

    // Face 2
    aFace.addNode(aNode5);
    aFace.addNode(aNode6);
    aFace.addNode(aNode7);
    aFace.addNode(aNode8);
    aFace.close();

    aControlVolume.addFace(1, aFace);
    aFace.reset();

    // Face 3
    aFace.addNode(aNode5);
    aFace.addNode(aNode1);
    aFace.addNode(aNode2);
    aFace.addNode(aNode6);
    aFace.close();

    aControlVolume.addFace(2, aFace);
    aFace.reset();

    // Face 4
    aFace.addNode(aNode8);
    aFace.addNode(aNode7);
    aFace.addNode(aNode3);
    aFace.addNode(aNode4);
    aFace.close();

    aControlVolume.addFace(3, aFace);
    aFace.reset();

    // Face 5
    aFace.addNode(aNode6);
    aFace.addNode(aNode2);
    aFace.addNode(aNode3);
    aFace.addNode(aNode7);
    aFace.close();
    aControlVolume.addFace(4, aFace);
    aFace.reset();

    // Face 6
    aFace.addNode(aNode1);
    aFace.addNode(aNode5);
    aFace.addNode(aNode8);
    aFace.addNode(aNode4);
    aFace.close();

    aControlVolume.addFace(5, aFace);
    aFace.reset();

    aControlVolume.calculateVolume(true);

    EXPECT_NEAR(1.0, aControlVolume.volume(), 1E-12);

    CGeoNormal<decimal> aNormal(1.0, 0.0, 0.0);

    decimal volumeFractionAct;
    Integer nIterations;

    aControlVolume.setClippedFaceId(999);
    aControlVolume.clip(aNormal, 1E-9, volumeFractionAct, nIterations, 50, 1E-12);

    EXPECT_NEAR(1E-9, volumeFractionAct, 1E-10);

}

