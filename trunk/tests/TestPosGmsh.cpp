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
#include "PosGmsh.hpp"

using namespace ENigMA::pde;
using namespace ENigMA::post;

class CTestPosGmsh : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestPosGmsh, line)
{

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<Decimal> aLine(aPoint1, aPoint2);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;

    aBasicMesher.generate(aLine, nu);

    EXPECT_EQ(nu, aBasicMesher.mesh().nbElements());

    CPdeField<Decimal> T1, T2;
    CPosGmsh<Decimal> aPosGmsh1, aPosGmsh2;

    T1.setMesh(aBasicMesher.mesh());
    aPosGmsh1.save(T1, "line.msh", "line");
    aPosGmsh2.load(T2, "line.msh");

    EXPECT_EQ(T1.mesh().nbNodes(), T2.mesh().nbNodes());
    EXPECT_EQ(T1.mesh().nbElements(), T2.mesh().nbElements());

}

TEST_F(CTestPosGmsh, rectangle1)
{

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 0.5, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 0.5, 0.0);

    CGeoQuadrilateral<Decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;

    aBasicMesher.generate(aQuadrilateral, nu, nv);

    EXPECT_EQ(nu * nv, aBasicMesher.mesh().nbElements());

    CPdeField<Decimal> T1, T2;
    CPosGmsh<Decimal> aPosGmsh1, aPosGmsh2;

    T1.setMesh(aBasicMesher.mesh());
    aPosGmsh1.save(T1, "rectangle1.msh", "rectangle1");
    aPosGmsh2.load(T2, "rectangle1.msh");

    EXPECT_EQ(T1.mesh().nbNodes(), T2.mesh().nbNodes());
    EXPECT_EQ(T1.mesh().nbElements(), T2.mesh().nbElements());

}

TEST_F(CTestPosGmsh, rectangle2)
{

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 0.5, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 0.5, 0.0);

    CGeoQuadrilateral<Decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;

    aBasicMesher.generate(aQuadrilateral, nu, nv, true);

    EXPECT_EQ(nu * nv * 2, aBasicMesher.mesh().nbElements());

    CPdeField<Decimal> T1, T2;
    CPosGmsh<Decimal> aPosGmsh1, aPosGmsh2;

    T1.setMesh(aBasicMesher.mesh());
    aPosGmsh1.save(T1, "rectangle2.msh", "rectangle2");
    aPosGmsh2.load(T2, "rectangle2.msh");

    EXPECT_EQ(T1.mesh().nbNodes(), T2.mesh().nbNodes());
    EXPECT_EQ(T1.mesh().nbElements(), T2.mesh().nbElements());

}

TEST_F(CTestPosGmsh, box1)
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

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    EXPECT_EQ(nu * nv * nw, aBasicMesher.mesh().nbElements());

    CPdeField<Decimal> T1, T2;
    CPosGmsh<Decimal> aPosGmsh1, aPosGmsh2;

    T1.setMesh(aBasicMesher.mesh());
    aPosGmsh1.save(T1, "box1.msh", "box1");
    aPosGmsh2.load(T2, "box1.msh");

    EXPECT_EQ(T1.mesh().nbNodes(), T2.mesh().nbNodes());
    EXPECT_EQ(T1.mesh().nbElements(), T2.mesh().nbElements());

}

TEST_F(CTestPosGmsh, box2)
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

    CPdeField<Decimal> T1, T2;
    CPosGmsh<Decimal> aPosGmsh1, aPosGmsh2;

    T1.setMesh(aBasicMesher.mesh());
    aPosGmsh1.save(T1, "box2.msh", "box2");
    aPosGmsh2.load(T2, "box2.msh");

    EXPECT_EQ(T1.mesh().nbNodes(), T2.mesh().nbNodes());
    EXPECT_EQ(T1.mesh().nbElements(), T2.mesh().nbElements());

}
