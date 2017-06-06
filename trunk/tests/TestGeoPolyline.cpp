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

#include "GeoPolyline.hpp"

using namespace ENigMA::geometry;

class CTestGeoPolyline : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoPolyline, create) {

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aPoint4(0.0, 1.0, 0.0);

    CGeoLine<decimal> aLine1(aPoint1, aPoint2);
    CGeoLine<decimal> aLine2(aPoint2, aPoint3);
    CGeoLine<decimal> aLine3(aPoint3, aPoint4);
    CGeoLine<decimal> aLine4(aPoint4, aPoint1);

    CGeoLineList<decimal> aLineList;

    aLineList.addLine(aLine1);
    aLineList.addLine(aLine2);
    aLineList.addLine(aLine3);
    aLineList.addLine(aLine4);

    CGeoPolyline<decimal> aPolyline(aLineList);

    EXPECT_EQ(4, aPolyline.nbLines());
    EXPECT_EQ(5, aPolyline.nbVertices());

}

TEST_F(CTestGeoPolyline, length) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);

    aPolyline.calculateLength();

    EXPECT_EQ(2.0, aPolyline.length());

}

TEST_F(CTestGeoPolyline, isClosed) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);

    aPolyline.close();

    EXPECT_TRUE(aPolyline.isClosed());

}

TEST_F(CTestGeoPolyline, nbLines) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);

    EXPECT_EQ(3, aPolyline.nbLines());

}

