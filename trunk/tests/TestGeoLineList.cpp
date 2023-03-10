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

#include "GeoLineList.hpp"

using namespace ENigMA::geometry;

class CTestGeoLineList : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoLineList, sort) {

    CGeoCoordinate<Decimal> aPoint1( 0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2( 1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint3( 3.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint4( 7.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint5(15.0, 0.0, 0.0);

    CGeoLine<Decimal> aLine1(aPoint1, aPoint2);
    CGeoLine<Decimal> aLine2(aPoint4, aPoint5);
    CGeoLine<Decimal> aLine3(aPoint3, aPoint4);
    CGeoLine<Decimal> aLine4(aPoint2, aPoint3);

    CGeoLineList<Decimal> aLineList;

    aLineList.addLine(aLine1);
    aLineList.addLine(aLine2);
    aLineList.addLine(aLine3);
    aLineList.addLine(aLine4);

    aLineList.sort(1E-12);

    aLineList.calculateLength(true);

    EXPECT_EQ(1.0, aLineList.line(0).length());
    EXPECT_EQ(2.0, aLineList.line(1).length());
    EXPECT_EQ(4.0, aLineList.line(2).length());
    EXPECT_EQ(8.0, aLineList.line(3).length());

}

