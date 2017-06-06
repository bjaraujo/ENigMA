// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#include "gtest/gtest.h"

#include "TypeDef.hpp"

#include "MshBasicMesher.hpp"
#include "MshExtrudedMesher.hpp"
#include "PdeField.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::pde;

class CTestMshExtrudedMesher : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestMshExtrudedMesher, extrudeTriangles) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoQuadrilateral<decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 4;

    aBasicMesher.generate(aQuadrilateral, nu, nv, true);

    CMshExtrudedMesher<decimal> anExtrudedMesher;

    CMshMesh<decimal>& aPlanarMesh = aBasicMesher.mesh();

    const decimal dw = 0.25;

    anExtrudedMesher.generate(aPlanarMesh, nw, dw, 1E-6);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anExtrudedMesher.mesh());
    aPosGmsh.save(T, "extruded_tris.msh", "prisms");

    EXPECT_EQ(nu * nv * nw * 2, anExtrudedMesher.mesh().nbElements());

}

TEST_F(CTestMshExtrudedMesher, extrudeQuadrilaterals) {

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoQuadrilateral<decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 4;

    aBasicMesher.generate(aQuadrilateral, nu, nv, false);

    CMshExtrudedMesher<decimal> anExtrudedMesher;

    CMshMesh<decimal>& aPlanarMesh = aBasicMesher.mesh();

    const decimal dw = 0.25;

    anExtrudedMesher.generate(aPlanarMesh, nw, dw, 1E-6);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anExtrudedMesher.mesh());
    aPosGmsh.save(T, "extruded_quads.msh", "hex");

    EXPECT_EQ(nu * nv * nw, anExtrudedMesher.mesh().nbElements());

}
