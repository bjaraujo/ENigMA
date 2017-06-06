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

#include "GeoHashGrid.hpp"

using namespace ENigMA::geometry;

class TestGeoHashGrid : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(TestGeoHashGrid, find1) {

    CGeoHashGrid<decimal> aHashGrid;
    
    Integer anIdToFind = 33;
    decimal x = 0.5;
    decimal y = 0.5;
    decimal z = 0.5;

    CGeoCoordinate<decimal> aCoord1(0.4, 0.5, 0.4);  // 0
    CGeoCoordinate<decimal> aCoord2(0.5, 0.6, 0.5);  // 1
    CGeoCoordinate<decimal> aCoord3(0.8, 0.1, 0.5);  // 2
    CGeoCoordinate<decimal> aCoord4(0.8, 0.1, 0.5);  // 3
    CGeoCoordinate<decimal> aCoord5(0.8, 0.1, 0.5);  // 4

    aHashGrid.addGeometricObject(12, aCoord1);
    aHashGrid.addGeometricObject(14, aCoord2);
    aHashGrid.addGeometricObject(15, aCoord3);
    aHashGrid.addGeometricObject(67, aCoord4);
    aHashGrid.addGeometricObject(89, aCoord5);

    for (int i = 0; i < 20; ++i)
    {
        CGeoCoordinate<decimal> aCoord(x + i * 0.1, y - i * 0.1, z + i * 0.1);
        aHashGrid.addGeometricObject(100 + i, aCoord);
    }

    CGeoCoordinate<decimal> aCoordToFind(x, y, z);
    aHashGrid.addGeometricObject(anIdToFind, aCoordToFind);

    for (int i = 20; i < 1000; ++i)
    {
        CGeoCoordinate<decimal> aCoord(x + (i - 20) * 0.1, y - (i - 20) * 0.1, z + (i - 20) * 0.1);
        aHashGrid.addGeometricObject(100 + i, aCoord);
    }

    aHashGrid.build();

    std::vector<Integer> sCoords;

    CGeoCoordinate<decimal> aCoord(x, y, z);
    aHashGrid.find(sCoords, aCoord, 1E-6);

    EXPECT_NE(sCoords.end(), std::find(sCoords.begin(), sCoords.end(), anIdToFind));

}

TEST_F(TestGeoHashGrid, find2) {

    CGeoHashGrid<decimal> aHashGrid;
    
    decimal d = 1.0;

    decimal x = 0.8*d;
    decimal y = 0.1*d;
    decimal z = 0.0*d;

    CGeoCoordinate<decimal> aCoord1(0.4*d, 0.5*d, 0.0);  // 0
    CGeoCoordinate<decimal> aCoord2(0.5*d, 0.6*d, 0.0);  // 1
    CGeoCoordinate<decimal> aCoord3(0.8*d, 0.1*d, 0.0);  // 2
    CGeoCoordinate<decimal> aCoord4(0.8*d, 0.1*d, 0.0);  // 3
    CGeoCoordinate<decimal> aCoord5(0.8*d, 0.1*d, 0.0);  // 4
    CGeoCoordinate<decimal> aCoord6(0.8*d, 0.2*d, 0.0);  // 5

    aHashGrid.addGeometricObject(12, aCoord1);
    aHashGrid.addGeometricObject(14, aCoord2);
    aHashGrid.addGeometricObject(15, aCoord3);
    aHashGrid.addGeometricObject(67, aCoord4);
    aHashGrid.addGeometricObject(89, aCoord5);
    aHashGrid.addGeometricObject(96, aCoord6);

    aHashGrid.build();

    std::vector<Integer> sCoords;

    CGeoCoordinate<decimal> aCoord(x, y, z);
    aHashGrid.find(sCoords, aCoord, 1E-6);

    EXPECT_EQ(3, sCoords.size());

}
