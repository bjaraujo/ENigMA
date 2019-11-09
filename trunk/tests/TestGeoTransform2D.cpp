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

#include "GeoCoordinate.hpp"
#include "GeoVector.hpp"
#include "GeoTriangle.hpp"

using namespace ENigMA::geometry;

class CTestGeoTransform2D : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoTransform2D, project2D) {

    CGeoCoordinate<decimal> aVertex1(+2.0, -3.0, +4.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +0.2, -3.0);
    CGeoCoordinate<decimal> aVertex3(-1.1, -1.0, -1.2);

    CGeoTriangle<decimal> aTriangle;

    aTriangle.addVertex(aVertex1);
    aTriangle.addVertex(aVertex2);
    aTriangle.addVertex(aVertex3);

    aTriangle.calculateArea();

    decimal a1 = aTriangle.area();

    CGeoVector<decimal> aVector1 = aVertex2 - aVertex1;
    CGeoVector<decimal> aVector2 = aVertex3 - aVertex1;

    CGeoVector<decimal> vx = aVector1;
    vx.normalize();

    CGeoVector<decimal> vz = vx.cross(aVector2);
    vz.normalize();

    CGeoVector<decimal> vy = vz.cross(vx);
    vy.normalize();

    CGeoCoordinateSystem<decimal> aCoordinateSystem;
    aCoordinateSystem.col(0) << vx;
    aCoordinateSystem.col(1) << vy;
    aCoordinateSystem.col(2) << vz;

    aCoordinateSystem.transposeInPlace();

    aVertex1.transform(aCoordinateSystem);
    aVertex2.transform(aCoordinateSystem);
    aVertex3.transform(aCoordinateSystem);

    /*
    std::cout << aVertex1 << std::endl;
    std::cout << aVertex2 << std::endl;
    std::cout << aVertex3 << std::endl;
    */

    aVertex2.z() -= aVertex1.z();
    aVertex3.z() -= aVertex1.z();
    aVertex1.z() = 0.0;

    aTriangle.calculateArea();

    decimal a2 = aTriangle.area();

    EXPECT_NEAR(0.0, aVertex1.z(), 1E-6);
    EXPECT_NEAR(0.0, aVertex2.z(), 1E-6);
    EXPECT_NEAR(0.0, aVertex3.z(), 1E-6);

    EXPECT_EQ(a1, a2);

}
