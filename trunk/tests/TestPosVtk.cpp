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
#include "PosVtk.hpp"

using namespace ENigMA::pde;
using namespace ENigMA::post;

class CTestPosVtk : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestPosVtk, box)
{

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 0.5, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 0.5, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, 0.5);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, 0.5);
    CGeoCoordinate<Decimal> aVertex7(1.0, 0.5, 0.5);
    CGeoCoordinate<Decimal> aVertex8(0.0, 0.5, 0.5);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 2;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    EXPECT_EQ(nu * nv * nw * 6, aBasicMesher.mesh().nbElements());

    CPdeField<Decimal> T;
    CPosVtk<Decimal> aPosVtk;

    T.setMesh(aBasicMesher.mesh());
    aPosVtk.save(T, "box.vtk");

    EXPECT_EQ(T.mesh().nbNodes(), 36);

}
