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

#include "AnaTemperature.hpp"

using namespace ENigMA::analytical;

class CTestAnaTemperature : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestAnaTemperature, steadyStateHeatConduction1D)
{

    CAnaTemperature<decimal> aAnaTemperature;

    decimal x, T;

    x = 0.5;
    aAnaTemperature.steadyStateHeatConduction1D(x, T);
    EXPECT_EQ(x, T);

}

TEST_F(CTestAnaTemperature, steadyStateHeatConduction2D)
{

    CAnaTemperature<decimal> aAnaTemperature;

    decimal x, y, T;

    x = 0.5;
    y = 0.5;
    aAnaTemperature.steadyStateHeatConduction2D(x, y, T);
    EXPECT_NEAR(0.25, T, 1E-2);

}

TEST_F(CTestAnaTemperature, steadyStateHeatConduction3D)
{

    CAnaTemperature<decimal> aAnaTemperature;

    decimal x, y, z, T;

    x = 0.5;
    y = 0.5;
    z = 0.5;
    aAnaTemperature.steadyStateHeatConduction3D(x, y, z, T);
    EXPECT_NEAR(0.25, T, 1E-2);

    x = 0.0;
    y = 0.5;
    z = 0.5;
    aAnaTemperature.steadyStateHeatConduction3D(x, y, z, T);
    EXPECT_NEAR(0.5, T, 1E-1);

    x = 1.0;
    y = 0.5;
    z = 0.5;
    aAnaTemperature.steadyStateHeatConduction3D(x, y, z, T);
    EXPECT_NEAR(1.0, T, 1E-1);

}

/*
TODO:

TEST_F(CTestAnaTemperature, steadyStateHeatConvectionRadiation1D)
{

    CAnaTemperature<decimal> aAnaTemperature;

    decimal x, Tb, h, e, k, perimeter, sectionArea, T;
    
    x = 0.0;
    Tb = 1.0;
    h = 1E-12;
    e = 3.0/2.0/5.6704E-8;
    k = 1.0;
    perimeter = 1.0;
    sectionArea = 1.0;

    aAnaTemperature.steadyStateHeatConvectionRadiation1D(x, Tb, h, e, k, perimeter, sectionArea, T);

    EXPECT_EQ(1, 2);

}
*/
