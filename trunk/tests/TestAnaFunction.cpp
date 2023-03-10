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

#include "AnaFunction.hpp"

using namespace ENigMA::analytical;

class CTestAnaFunction : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestAnaFunction, evaluate)
{

    Decimal x;
    Decimal f;

    CAnaFunction<Decimal> aAnaFunction;

    x = 1.0;
    f = 1.0;
    aAnaFunction.set("f * cos(x)");
    aAnaFunction.defineVariable("x", x);
    aAnaFunction.defineVariable("f", f);
    EXPECT_NEAR(f * cos(x), aAnaFunction.evaluate(), 1E-6);

    x = 0.7;
    aAnaFunction.set("cos(x)*sin(x)/2");
    EXPECT_NEAR(cos(x)*sin(x)/2, aAnaFunction.evaluate(), 1E-6);

    Integer nIterBisect, nIterBrent;

    aAnaFunction.set("2*x-6");
    EXPECT_NEAR(3.0, aAnaFunction.bisection("x", -10.0, 10.0, nIterBisect, 100, 1E-12), 1E-6);

    aAnaFunction.set("2*x-6");
    EXPECT_NEAR(3.0, aAnaFunction.brent("x", -10.0, 10.0, nIterBrent, 100, 1E-12), 1E-6);

    EXPECT_LT(nIterBrent, nIterBisect);

    aAnaFunction.set("x^3-x-2");
    EXPECT_NEAR(1.5213623, aAnaFunction.bisection("x", -10.0, 10.0, nIterBisect, 100, 1E-12), 1E-3);

    aAnaFunction.set("x^3-x-2");
    EXPECT_NEAR(1.5213623, aAnaFunction.brent("x", -10.0, 10.0, nIterBrent, 100, 1E-12), 1E-3);

    EXPECT_LT(nIterBrent, nIterBisect);

    Integer nIterRoot;

    aAnaFunction.set("x^2-6*x+9");
    EXPECT_NEAR(3.0, aAnaFunction.root("x", 2.0, 4.0, nIterRoot, 100, 1E-12), 1E-3);

}

