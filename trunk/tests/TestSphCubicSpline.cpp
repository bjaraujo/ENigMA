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
#include "GeoVertexList.hpp"
#include "SphCubicSpline.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::sph;

class CTestSphCubicSpline : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestSphCubicSpline, gradient1)
{

    CGeoCoordinate<decimal> aVertex1(+0.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex3(-1.0, +0.0, +0.0);

    CGeoVertexList<decimal> aVertexList;

    aVertexList.addVertex(aVertex1);
    aVertexList.addVertex(aVertex2);
    aVertexList.addVertex(aVertex3);

    std::vector<decimal> sValues;

    sValues.emplace_back(0);
    sValues.emplace_back(1);
    sValues.emplace_back(1);

    CSphCubicSpline<decimal> aKernel(1);

    decimal gradx = 0.0;

    for (Integer i = 1; i < aVertexList.nbVertices(); i++)
    {

        CGeoVector<decimal> r = aVertexList.vertex(i) - aVertexList.vertex(0);

        gradx += sValues[i] * aKernel.gradientW(r, 1.0).x();

    }

    EXPECT_NEAR(0.0, gradx, 1E-12);

    sValues[1] = +1;
    sValues[2] = -1;

    gradx = 0.0;

    for (Integer i = 1; i < aVertexList.nbVertices(); i++)
    {

        CGeoVector<decimal> r = aVertexList.vertex(i) - aVertexList.vertex(0);

        gradx += sValues[i] * aKernel.gradientW(r, 1.0).x();

    }

    EXPECT_NEAR(1.0, gradx, 0.5);

    sValues[1] = +3;
    sValues[2] = -3;

    gradx = 0.0;

    for (Integer i = 1; i < aVertexList.nbVertices(); i++)
    {

        CGeoVector<decimal> r = aVertexList.vertex(i) - aVertexList.vertex(0);

        gradx += sValues[i] * aKernel.gradientW(r, 1.0).x();

    }

    EXPECT_NEAR(3.0, gradx, 1.5);

}

TEST_F(CTestSphCubicSpline, gradient2)
{

    CGeoCoordinate<decimal> aVertex1(+0.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex3(-1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex4(+0.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex5(+0.0, -1.0, +0.0);
    CGeoCoordinate<decimal> aVertex6(+0.5, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex7(-0.5, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex8(+0.0, +0.5, +0.0);
    CGeoCoordinate<decimal> aVertex9(+0.0, -0.5, +0.0);

    CGeoVertexList<decimal> aVertexList;

    aVertexList.addVertex(aVertex1);
    aVertexList.addVertex(aVertex2);
    aVertexList.addVertex(aVertex3);
    aVertexList.addVertex(aVertex4);
    aVertexList.addVertex(aVertex5);
    aVertexList.addVertex(aVertex6);
    aVertexList.addVertex(aVertex7);
    aVertexList.addVertex(aVertex8);
    aVertexList.addVertex(aVertex9);

    std::vector<decimal> sValues;

    sValues.emplace_back(0);
    sValues.emplace_back(1);
    sValues.emplace_back(1);
    sValues.emplace_back(1);
    sValues.emplace_back(1);
    sValues.emplace_back(1);
    sValues.emplace_back(1);
    sValues.emplace_back(1);
    sValues.emplace_back(1);

    CSphCubicSpline<decimal> aKernel(2);

    decimal gradx = 0.0;
    decimal grady = 0.0;

    for (Integer i = 1; i < aVertexList.nbVertices(); i++)
    {

        CGeoVector<decimal> r = aVertexList.vertex(i) - aVertexList.vertex(0);

        gradx += sValues[i] * aKernel.gradientW(r, 1.0).x();
        grady += sValues[i] * aKernel.gradientW(r, 1.0).y();

    }

    EXPECT_NEAR(0.0, gradx, 1E-12);
    EXPECT_NEAR(0.0, grady, 1E-12);

    sValues[1] = +1;
    sValues[2] = -1;
    sValues[3] = +1;
    sValues[4] = -1;
    sValues[5] = +0.5;
    sValues[6] = -0.5;
    sValues[7] = +0.5;
    sValues[8] = -0.5;

    gradx = 0.0;
    grady = 0.0;

    for (Integer i = 1; i < aVertexList.nbVertices(); i++)
    {

        CGeoVector<decimal> r = aVertexList.vertex(i) - aVertexList.vertex(0);

        gradx += sValues[i] * aKernel.gradientW(r, 1.0).x();
        grady += sValues[i] * aKernel.gradientW(r, 1.0).y();

    }

    EXPECT_NEAR(1.0, gradx, 0.5);
    EXPECT_NEAR(1.0, grady, 0.5);

    sValues[1] = +3;
    sValues[2] = -3;
    sValues[3] = +3;
    sValues[4] = -3;
    sValues[5] = +1.5;
    sValues[6] = -1.5;
    sValues[7] = +1.5;
    sValues[8] = -1.5;

    gradx = 0.0;
    grady = 0.0;

    for (Integer i = 1; i < aVertexList.nbVertices(); i++)
    {

        CGeoVector<decimal> r = aVertexList.vertex(i) - aVertexList.vertex(0);

        gradx += sValues[i] * aKernel.gradientW(r, 1.0).x();
        grady += sValues[i] * aKernel.gradientW(r, 1.0).y();

    }

    EXPECT_NEAR(3.0, gradx, 1.5);
    EXPECT_NEAR(3.0, grady, 1.5);

}

TEST_F(CTestSphCubicSpline, laplacian)
{

    CGeoCoordinate<decimal> aVertex1(+0.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex2(+1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex3(-1.0, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex4(+0.0, +1.0, +0.0);
    CGeoCoordinate<decimal> aVertex5(+0.0, -1.0, +0.0);
    CGeoCoordinate<decimal> aVertex6(+0.5, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex7(-0.5, +0.0, +0.0);
    CGeoCoordinate<decimal> aVertex8(+0.0, +0.5, +0.0);
    CGeoCoordinate<decimal> aVertex9(+0.0, -0.5, +0.0);

    CGeoVertexList<decimal> aVertexList;

    aVertexList.addVertex(aVertex1);
    aVertexList.addVertex(aVertex2);
    aVertexList.addVertex(aVertex3);
    aVertexList.addVertex(aVertex4);
    aVertexList.addVertex(aVertex5);
    aVertexList.addVertex(aVertex6);
    aVertexList.addVertex(aVertex7);
    aVertexList.addVertex(aVertex8);
    aVertexList.addVertex(aVertex9);

    std::vector<decimal> sValues;

    sValues.emplace_back(0);
    sValues.emplace_back(+1);
    sValues.emplace_back(-1);
    sValues.emplace_back(+1);
    sValues.emplace_back(-1);
    sValues.emplace_back(+0.5);
    sValues.emplace_back(-0.5);
    sValues.emplace_back(+0.5);
    sValues.emplace_back(-0.5);

    CSphCubicSpline<decimal> aKernel(2);

    decimal lap = 0.0;

    for (Integer i = 1; i < aVertexList.nbVertices(); i++)
    {

        CGeoVector<decimal> r = aVertexList.vertex(i) - aVertexList.vertex(0);

        lap += (sValues[i] - sValues[0]) * aKernel.laplacianW(r, 2.0);

    }

    EXPECT_NEAR(0.0, lap, 1E-12);

    sValues[0] = +0;
    sValues[1] = +3;
    sValues[2] = -1;
    sValues[3] = +3;
    sValues[4] = -1;
    sValues[5] = +0.5;
    sValues[6] = -0.5;
    sValues[7] = +0.5;
    sValues[8] = -0.5;

    lap = 0.0;

    for (Integer i = 1; i < aVertexList.nbVertices(); i++)
    {

        CGeoVector<decimal> r = aVertexList.vertex(i) - aVertexList.vertex(0);

        lap += (sValues[i] - sValues[0]) * aKernel.laplacianW(r, 2.0);

    }

    EXPECT_NEAR(1.5, lap, 0.5);

}
