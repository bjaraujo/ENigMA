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

#include "GeoTriangularPrism.hpp"
#include "GeoHexahedron.hpp"
#include "GeoPolyhedron.hpp"
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

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoQuadrilateral<Decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 4;

    aBasicMesher.generate(aQuadrilateral, nu, nv, true);

    CMshExtrudedMesher<Decimal> anExtrudedMesher;

    CMshMesh<Decimal>& aPlanarMesh = aBasicMesher.mesh();

    CPdeField<Decimal> T;
    CPosGmsh<Decimal> aPosGmsh;

    T.setMesh(aPlanarMesh);
    aPosGmsh.save(T, "planar_tris.msh", "tris");

    const Decimal dw = 0.25;

    for (Integer i = 0; i < nw; i++)
    {
        Decimal aDelta = dw / nw;
        anExtrudedMesher.generate(aPlanarMesh, aDelta, 1E-6);
    }

    CMshMesh<Decimal> aVolumeMesh;
    aVolumeMesh = anExtrudedMesher.mesh();

    T.setMesh(aVolumeMesh);
    aPosGmsh.save(T, "extruded_tris.msh", "prisms");

    EXPECT_EQ(nu * nv * nw * 2, aVolumeMesh.nbElements());

    CMshElement<Decimal> anElement = aVolumeMesh.element(0);

    CGeoTriangularPrism<Decimal> aPrism;
    aPrism.addVertex(aVolumeMesh.node(anElement.nodeId(0)));
    aPrism.addVertex(aVolumeMesh.node(anElement.nodeId(1)));
    aPrism.addVertex(aVolumeMesh.node(anElement.nodeId(2)));
    aPrism.addVertex(aVolumeMesh.node(anElement.nodeId(3)));
    aPrism.addVertex(aVolumeMesh.node(anElement.nodeId(4)));
    aPrism.addVertex(aVolumeMesh.node(anElement.nodeId(5)));
    aPrism.calculateVolume();

    EXPECT_GT(aPrism.volume(), 0.0);

    CGeoPolyhedron<Decimal> aPolyhedron(aPrism);
    aPolyhedron.calculateVolume();

    EXPECT_GT(aPolyhedron.volume(), 0.0);
}

TEST_F(CTestMshExtrudedMesher, extrudeQuadrilaterals) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoQuadrilateral<Decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 4;

    aBasicMesher.generate(aQuadrilateral, nu, nv, false);

    CMshExtrudedMesher<Decimal> anExtrudedMesher;

    CMshMesh<Decimal>& aPlanarMesh = aBasicMesher.mesh();

    const Decimal dw = 0.25;

    for (Integer i = 0; i < nw; i++)
    {
        Decimal aDelta = dw / nw;
        anExtrudedMesher.generate(aPlanarMesh, aDelta, 1E-6);
    }

    CMshMesh<Decimal> aVolumeMesh;
    aVolumeMesh = anExtrudedMesher.mesh();

    CPdeField<Decimal> T;
    CPosGmsh<Decimal> aPosGmsh;

    T.setMesh(aVolumeMesh);
    aPosGmsh.save(T, "extruded_quads.msh", "hex");

    EXPECT_EQ(nu * nv * nw, aVolumeMesh.nbElements());

    CMshElement<Decimal> anElement = aVolumeMesh.element(0);

    CGeoHexahedron<Decimal> aBrick;
    aBrick.addVertex(aVolumeMesh.node(anElement.nodeId(0)));
    aBrick.addVertex(aVolumeMesh.node(anElement.nodeId(1)));
    aBrick.addVertex(aVolumeMesh.node(anElement.nodeId(2)));
    aBrick.addVertex(aVolumeMesh.node(anElement.nodeId(3)));
    aBrick.addVertex(aVolumeMesh.node(anElement.nodeId(4)));
    aBrick.addVertex(aVolumeMesh.node(anElement.nodeId(5)));
    aBrick.addVertex(aVolumeMesh.node(anElement.nodeId(6)));
    aBrick.addVertex(aVolumeMesh.node(anElement.nodeId(7)));
    aBrick.calculateVolume();

    EXPECT_GT(aBrick.volume(), 0.0);

    CGeoPolyhedron<Decimal> aPolyhedron(aBrick);
    aPolyhedron.calculateVolume();

    EXPECT_GT(aPolyhedron.volume(), 0.0);
}
