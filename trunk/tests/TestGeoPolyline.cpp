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

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4(0.0, 1.0, 0.0);

    CGeoLine<Decimal> aLine1(aPoint1, aPoint2);
    CGeoLine<Decimal> aLine2(aPoint2, aPoint3);
    CGeoLine<Decimal> aLine3(aPoint3, aPoint4);
    CGeoLine<Decimal> aLine4(aPoint4, aPoint1);

    CGeoLineList<Decimal> aLineList;

    aLineList.addLine(aLine1);
    aLineList.addLine(aLine2);
    aLineList.addLine(aLine3);
    aLineList.addLine(aLine4);

    CGeoPolyline<Decimal> aPolyline(aLineList);

    EXPECT_EQ(4, aPolyline.nbLines());
    EXPECT_EQ(5, aPolyline.nbVertices());

}

TEST_F(CTestGeoPolyline, length) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);

    CGeoPolyline<Decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);

    aPolyline.calculateLength();

    EXPECT_EQ(2.0, aPolyline.length());

}

TEST_F(CTestGeoPolyline, isClosed) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<Decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);

    aPolyline.close();

    EXPECT_TRUE(aPolyline.isClosed());

}

TEST_F(CTestGeoPolyline, nbLines) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoPolyline<Decimal> aPolyline;

    aPolyline.addVertex(aVertex1);
    aPolyline.addVertex(aVertex2);
    aPolyline.addVertex(aVertex3);
    aPolyline.addVertex(aVertex4);

    EXPECT_EQ(3, aPolyline.nbLines());

}

