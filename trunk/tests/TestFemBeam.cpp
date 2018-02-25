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

#include "FemBeam.hpp"

using namespace ENigMA::fem;

class CTestFemBeam : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFemBeam, update) {

    CFemBeam<decimal, 2, 1, 1> aBeam;

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 1.0, 1.0);

    aBeam.setStartPoint(aPoint1);
    aBeam.setEndPoint(aPoint2);

    aBeam.calculateLength();

    aBeam.setTransient(true);

    aBeam.setDt(0.1);

    aBeam.update();

    //std::cout << aBeam.laplacian << std::endl;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {

            decimal s;

            if (i == j) s = -1; else s = +1;

            EXPECT_NEAR(1.0 / sqrt(3.0) * s, aBeam.laplacian(i, j), 1E-6);

        }
    }

    //std::cout << aBeam.ddt << std::endl;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {

            decimal p;

            if (i == j) p = 1.0/3.0; else p = 1.0/6.0;

            EXPECT_NEAR(p * aBeam.length() * aBeam.sectionArea() / aBeam.dt(), aBeam.ddt(i, j), 1E-6);

        }
    }

}

