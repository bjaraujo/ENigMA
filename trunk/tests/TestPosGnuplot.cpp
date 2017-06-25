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

#include "MshBasicMesher.hpp"
#include "PdeField.hpp"
#include "PdeEquation.hpp"
#include "PosGnuplot.hpp"

using namespace ENigMA::pde;
using namespace ENigMA::post;

class CTestPosGnuplot : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestPosGnuplot, plot)
{

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<decimal> aLine(aPoint1, aPoint2);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 3;

    aBasicMesher.generate(aLine, nu);

    EXPECT_EQ(nu, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;
    CPosGnuplot<decimal> aPosGnuplot;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretLocation(DL_NODE);
    T.u.resize(T.mesh().nbNodes());

    aPosGnuplot.save(T, "plot.dat");

    EXPECT_EQ(T.u.size(), 4);

}

