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

#include "GeoAdtree.hpp"

using namespace ENigMA::geometry;

class TestGeoAdtree : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(TestGeoAdtree, find1) {

    CGeoBoundingBox<decimal> aBoundingBox1;
    CGeoBoundingBox<decimal> aBoundingBox2;
    CGeoBoundingBox<decimal> aBoundingBox3;

    CGeoCoordinate<decimal> aCoordinate1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate2(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate3(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate4(3.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate5(4.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate6(5.0, 1.0, 0.0);

    aBoundingBox1.addCoordinate(aCoordinate1);
    aBoundingBox1.addCoordinate(aCoordinate2);

    aBoundingBox2.addCoordinate(aCoordinate3);
    aBoundingBox2.addCoordinate(aCoordinate4);

    aBoundingBox3.addCoordinate(aCoordinate5);
    aBoundingBox3.addCoordinate(aCoordinate6);

    CGeoBoundingBox<decimal> aBoundingBox;

    aBoundingBox.reset();
    aBoundingBox.addCoordinate(aCoordinate1);
    aBoundingBox.addCoordinate(aCoordinate2);
    aBoundingBox.addCoordinate(aCoordinate3);
    aBoundingBox.addCoordinate(aCoordinate4);
    aBoundingBox.addCoordinate(aCoordinate5);
    aBoundingBox.addCoordinate(aCoordinate6);

    CGeoAdtree<decimal> anAdtree;

    anAdtree.set(aBoundingBox);
    anAdtree.addGeometricObject(12, aBoundingBox1);
    anAdtree.addGeometricObject(14, aBoundingBox2);
    anAdtree.addGeometricObject(15, aBoundingBox3);

    anAdtree.build(); // Optional

    CGeoCoordinate<decimal> aCoordinate7(2.5, -0.5, 0.0);
    CGeoCoordinate<decimal> aCoordinate8(4.5, +0.5, 0.0);

    aBoundingBox.reset();
    aBoundingBox.addCoordinate(aCoordinate7);
    aBoundingBox.addCoordinate(aCoordinate8);

    std::vector<Integer> sBoundingBoxes;
    anAdtree.find(sBoundingBoxes, aBoundingBox);

    EXPECT_EQ(2, sBoundingBoxes.size());

    EXPECT_NE(sBoundingBoxes.end(), std::find(sBoundingBoxes.begin(), sBoundingBoxes.end(), 14));
    EXPECT_NE(sBoundingBoxes.end(), std::find(sBoundingBoxes.begin(), sBoundingBoxes.end(), 15));

}

TEST_F(TestGeoAdtree, find2) {

    CGeoBoundingBox<decimal> aBoundingBox1;
    CGeoBoundingBox<decimal> aBoundingBox2;
    CGeoBoundingBox<decimal> aBoundingBox3;
    CGeoBoundingBox<decimal> aBoundingBox4;
    CGeoBoundingBox<decimal> aBoundingBox5;
    CGeoBoundingBox<decimal> aBoundingBox6;

    CGeoCoordinate<decimal> aCoordinate01(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate02(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate03(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate04(3.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate05(4.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate06(5.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate07(0.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate08(4.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate09(5.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate10(6.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate11(7.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate12(8.0, 1.0, 1.0);

    aBoundingBox1.addCoordinate(aCoordinate01);
    aBoundingBox1.addCoordinate(aCoordinate02);

    aBoundingBox2.addCoordinate(aCoordinate03);
    aBoundingBox2.addCoordinate(aCoordinate04);

    aBoundingBox3.addCoordinate(aCoordinate05);
    aBoundingBox3.addCoordinate(aCoordinate06);

    aBoundingBox4.addCoordinate(aCoordinate07);
    aBoundingBox4.addCoordinate(aCoordinate08);

    aBoundingBox5.addCoordinate(aCoordinate09);
    aBoundingBox5.addCoordinate(aCoordinate10);

    aBoundingBox6.addCoordinate(aCoordinate11);
    aBoundingBox6.addCoordinate(aCoordinate12);

    CGeoBoundingBox<decimal> aBoundingBox;

    aBoundingBox.reset();
    aBoundingBox.addCoordinate(aCoordinate01);
    aBoundingBox.addCoordinate(aCoordinate02);
    aBoundingBox.addCoordinate(aCoordinate03);
    aBoundingBox.addCoordinate(aCoordinate04);
    aBoundingBox.addCoordinate(aCoordinate05);
    aBoundingBox.addCoordinate(aCoordinate06);
    aBoundingBox.addCoordinate(aCoordinate07);
    aBoundingBox.addCoordinate(aCoordinate08);
    aBoundingBox.addCoordinate(aCoordinate09);
    aBoundingBox.addCoordinate(aCoordinate10);
    aBoundingBox.addCoordinate(aCoordinate11);
    aBoundingBox.addCoordinate(aCoordinate12);

    CGeoAdtree<decimal> anAdtree;

    anAdtree.set(aBoundingBox);
    anAdtree.addGeometricObject(12, aBoundingBox1);
    anAdtree.addGeometricObject(14, aBoundingBox2);
    anAdtree.addGeometricObject(15, aBoundingBox3);
    anAdtree.addGeometricObject(21, aBoundingBox4);
    anAdtree.addGeometricObject(22, aBoundingBox5);
    anAdtree.addGeometricObject(30, aBoundingBox6);

    anAdtree.build(); // Optional

    CGeoCoordinate<decimal> aCoordinate13(2.5, +0.5, +0.5);
    CGeoCoordinate<decimal> aCoordinate14(4.5, -0.5, -0.5);

    aBoundingBox.reset();
    aBoundingBox.addCoordinate(aCoordinate13);
    aBoundingBox.addCoordinate(aCoordinate14);

    std::vector<Integer> sBoundingBoxes;
    anAdtree.find(sBoundingBoxes, aBoundingBox);

    EXPECT_EQ(2, sBoundingBoxes.size());

    EXPECT_NE(sBoundingBoxes.end(), std::find(sBoundingBoxes.begin(), sBoundingBoxes.end(), 14));
    EXPECT_NE(sBoundingBoxes.end(), std::find(sBoundingBoxes.begin(), sBoundingBoxes.end(), 15));

}

TEST_F(TestGeoAdtree, find3) {

    CGeoBoundingBox<decimal> aBoundingBox1;
    CGeoBoundingBox<decimal> aBoundingBox2;
    CGeoBoundingBox<decimal> aBoundingBox3;
    CGeoBoundingBox<decimal> aBoundingBox4;
    CGeoBoundingBox<decimal> aBoundingBox5;
    CGeoBoundingBox<decimal> aBoundingBox6;

    CGeoCoordinate<decimal> aCoordinate01(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate02(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate03(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate04(3.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate05(4.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate06(5.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aCoordinate07(0.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate08(4.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate09(5.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate10(6.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate11(7.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate12(8.0, 1.0, 1.0);

    aBoundingBox1.addCoordinate(aCoordinate01);
    aBoundingBox1.addCoordinate(aCoordinate02);

    aBoundingBox2.addCoordinate(aCoordinate03);
    aBoundingBox2.addCoordinate(aCoordinate04);

    aBoundingBox3.addCoordinate(aCoordinate05);
    aBoundingBox3.addCoordinate(aCoordinate06);

    aBoundingBox4.addCoordinate(aCoordinate07);
    aBoundingBox4.addCoordinate(aCoordinate08);

    aBoundingBox5.addCoordinate(aCoordinate09);
    aBoundingBox5.addCoordinate(aCoordinate10);

    aBoundingBox6.addCoordinate(aCoordinate11);
    aBoundingBox6.addCoordinate(aCoordinate12);

    CGeoBoundingBox<decimal> aBoundingBox;

    aBoundingBox.reset();
    aBoundingBox.addCoordinate(aCoordinate01);
    aBoundingBox.addCoordinate(aCoordinate02);
    aBoundingBox.addCoordinate(aCoordinate03);
    aBoundingBox.addCoordinate(aCoordinate04);
    aBoundingBox.addCoordinate(aCoordinate05);
    aBoundingBox.addCoordinate(aCoordinate06);
    aBoundingBox.addCoordinate(aCoordinate07);
    aBoundingBox.addCoordinate(aCoordinate08);
    aBoundingBox.addCoordinate(aCoordinate09);
    aBoundingBox.addCoordinate(aCoordinate10);
    aBoundingBox.addCoordinate(aCoordinate11);
    aBoundingBox.addCoordinate(aCoordinate12);

    CGeoAdtree<decimal> anAdtree;

    anAdtree.set(aBoundingBox);
    anAdtree.addGeometricObject(12, aBoundingBox1);
    anAdtree.addGeometricObject(14, aBoundingBox2);
    anAdtree.addGeometricObject(15, aBoundingBox3);
    anAdtree.addGeometricObject(21, aBoundingBox4);
    anAdtree.addGeometricObject(22, aBoundingBox5);
    anAdtree.addGeometricObject(30, aBoundingBox6);

    anAdtree.build(); // Optional

    CGeoCoordinate<decimal> aCoordinate13(0.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aCoordinate14(0.0, 1.0, 2.0);

    aBoundingBox.reset();
    aBoundingBox.addCoordinate(aCoordinate13);
    aBoundingBox.addCoordinate(aCoordinate14);

    std::vector<Integer> sBoundingBoxes;
    anAdtree.find(sBoundingBoxes, aBoundingBox);

    EXPECT_EQ(1, sBoundingBoxes.size());

    EXPECT_NE(sBoundingBoxes.end(), std::find(sBoundingBoxes.begin(), sBoundingBoxes.end(), 21));

}
