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

#include "GeoOctree.hpp"

using namespace ENigMA::geometry;

class TestGeoOctree : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(TestGeoOctree, find) {

    CGeoOctree<decimal> anOctree;
    
    Integer anIdToFind = 33;
    decimal x = 0.5;
    decimal y = 0.5;
    decimal z = 0.5;

    CGeoCoordinate<decimal> aCoord1(0.4, 0.5, 0.4);  // 0
    CGeoCoordinate<decimal> aCoord2(0.5, 0.6, 0.5);  // 1
    CGeoCoordinate<decimal> aCoord3(0.8, 0.1, 0.5);  // 3
    CGeoCoordinate<decimal> aCoord4(0.8, 0.1, 0.5);  // 4
    CGeoCoordinate<decimal> aCoord5(0.8, 0.1, 0.5);  // 5

    anOctree.addGeometricObject(12, aCoord1);
    anOctree.addGeometricObject(14, aCoord2);
    anOctree.addGeometricObject(15, aCoord3);
    anOctree.addGeometricObject(67, aCoord4);
    anOctree.addGeometricObject(89, aCoord5);

    for (int i = 0; i < 20; ++i)
    {
        CGeoCoordinate<decimal> aCoord(x + i * 0.1, y - i * 0.1, z + i * 0.1);
        anOctree.addGeometricObject(100 + i, aCoord);
    }

    CGeoCoordinate<decimal> aCoordToFind(x, y, z);
    anOctree.addGeometricObject(anIdToFind, aCoordToFind);

    for (int i = 20; i < 1000; ++i)
    {
        CGeoCoordinate<decimal> aCoord(x + (i - 20) * 0.1, y - (i - 20) * 0.1, z + (i - 20) * 0.1);
        anOctree.addGeometricObject(100 + i, aCoord);
    }

    anOctree.build();

    std::vector<Integer> sCoords;

    CGeoCoordinate<decimal> aCoord(x, y, z);
    anOctree.find(sCoords, aCoord, 1E-6);

    EXPECT_NE(sCoords.end(), std::find(sCoords.begin(), sCoords.end(), anIdToFind));

}

