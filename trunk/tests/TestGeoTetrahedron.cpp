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

#include "GeoTetrahedron.hpp"

using namespace ENigMA::geometry;

class CTestGeoTetrahedron : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoTetrahedron, volume1) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.0, 1.0);

    CGeoTetrahedron<decimal> aTetrahedron;

    aTetrahedron.addVertex(aVertex1);
    aTetrahedron.addVertex(aVertex2);
    aTetrahedron.addVertex(aVertex3);
    aTetrahedron.addVertex(aVertex4);

    aTetrahedron.calculateVolume();

    EXPECT_NEAR(1.0/6.0, aTetrahedron.volume(), 1E-6);
}

TEST_F(CTestGeoTetrahedron, volume2) {

    CGeoCoordinate<decimal> aVertex1(0.0206083, -0.00442002, 0.0089);
    CGeoCoordinate<decimal> aVertex2(0.0214445, -0.004, 0.0079416);
    CGeoCoordinate<decimal> aVertex3(0.02147799300482443, -0.00452017265334567, 0.008806528084091029);
    CGeoCoordinate<decimal> aVertex4(0.0214445, -0.004, 0.0089);

    CGeoTetrahedron<decimal> aTetrahedron;

    aTetrahedron.addVertex(aVertex1);
    aTetrahedron.addVertex(aVertex2);
    aTetrahedron.addVertex(aVertex3);
    aTetrahedron.addVertex(aVertex4);

    aTetrahedron.calculateVolume();

    EXPECT_GT(aTetrahedron.volume(), 0.0);

}

TEST_F(CTestGeoTetrahedron, contains1) {

    CGeoCoordinate<decimal> aVertex1(1, 0.375, 0.25);
    CGeoCoordinate<decimal> aVertex2(1.125, 0.25, 0.375);
    CGeoCoordinate<decimal> aVertex3(1.08333, 0.208333, 0.125);
    CGeoCoordinate<decimal> aVertex4(1.08333, 0.125, 0.208333);

    CGeoTetrahedron<decimal> aTetrahedron;

    aTetrahedron.addVertex(aVertex1);
    aTetrahedron.addVertex(aVertex2);
    aTetrahedron.addVertex(aVertex3);
    aTetrahedron.addVertex(aVertex4);

    CGeoCoordinate<decimal> aPoint(1.09324, 0.218523, 0.209264);

    EXPECT_TRUE(aTetrahedron.contains(aPoint, 1E-8));

}

TEST_F(CTestGeoTetrahedron, contains2) {

    CGeoCoordinate<decimal> aVertex1(0.5,0.25,0);
    CGeoCoordinate<decimal> aVertex2(0.25,0.25,0.25);
    CGeoCoordinate<decimal> aVertex3(0.5,0.25,0.25);
    CGeoCoordinate<decimal> aVertex4(0.5,0,0);

    CGeoTetrahedron<decimal> aTetrahedron;

    aTetrahedron.addVertex(aVertex1);
    aTetrahedron.addVertex(aVertex2);
    aTetrahedron.addVertex(aVertex3);
    aTetrahedron.addVertex(aVertex4);

    CGeoCoordinate<decimal> aPoint(0.416667, 0.211077, 0.122256);

    EXPECT_TRUE(aTetrahedron.contains(aPoint, 1E-8));

}

TEST_F(CTestGeoTetrahedron, contains3) {

    CGeoCoordinate<decimal> aVertex1(0.75,0.25,0);
    CGeoCoordinate<decimal> aVertex2(0.5,0.25,0.25);
    CGeoCoordinate<decimal> aVertex3(0.75,0.25,0.25);
    CGeoCoordinate<decimal> aVertex4(0.75,0,0);

    CGeoTetrahedron<decimal> aTetrahedron;

    aTetrahedron.addVertex(aVertex1);
    aTetrahedron.addVertex(aVertex2);
    aTetrahedron.addVertex(aVertex3);
    aTetrahedron.addVertex(aVertex4);

    CGeoCoordinate<decimal> aPoint(0.666667, 0.211077, 0.122256);

    EXPECT_TRUE(aTetrahedron.contains(aPoint, 1E-8));

}
