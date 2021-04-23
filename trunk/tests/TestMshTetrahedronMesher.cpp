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
#include "MshTetrahedronMesher.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::pde;

class CTestMshTetrahedronMesher : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestMshTetrahedronMesher, mesh1) {

    const decimal d = 0.125;

    const Integer nu = 2;
    const Integer nv = 2;
    const Integer nw = 2;

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(nu * d, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(nu * d, nv * d, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, nv * d, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, nw * d);
    CGeoCoordinate<decimal> aVertex6(nu * d, 0.0, nw * d);
    CGeoCoordinate<decimal> aVertex7(nu * d, nv * d, nw * d);
    CGeoCoordinate<decimal> aVertex8(0.0, nv * d, nw * d);

    CGeoHexahedron<decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CMshBasicMesher<decimal> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    CMshMesh<decimal> aSurfaceMesh = aBasicMesher.mesh().extractBoundary(1E-6);

    aSurfaceMesh.generateFaces(1E-5);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tetra_surface1.msh", "tris");

    CMshTetrahedronMesher<decimal> aTetrahedronMesher;

    aTetrahedronMesher.generate(aSurfaceMesh, 999, d * 1.5, d * 0.1, d * 10.0, 1E-3);

    CMshMesh<decimal> aVolumeMesh;
    aVolumeMesh = aTetrahedronMesher.mesh();

    T.setMesh(aVolumeMesh);
    aPosGmsh.save(T, "tetra_volume1.msh", "tetras");

    EXPECT_EQ(52, aVolumeMesh.nbElements());

}

TEST_F(CTestMshTetrahedronMesher, mesh2) {

    const decimal d = 0.125;

    const Integer nu = 4;
    const Integer nv = 2;
    const Integer nw = 2;    

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(nu * d, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(nu * d, nv * d, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, nv * d, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, nw * d);
    CGeoCoordinate<decimal> aVertex6(nu * d, 0.0, nw * d);
    CGeoCoordinate<decimal> aVertex7(nu * d, nv * d, nw * d);
    CGeoCoordinate<decimal> aVertex8(0.0, nv * d, nw * d);

    CGeoHexahedron<decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CMshBasicMesher<decimal> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    CMshMesh<decimal> aSurfaceMesh = aBasicMesher.mesh().extractBoundary(1E-6);

    aSurfaceMesh.generateFaces(1E-5);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tetra_surface2.msh", "tris");

    CMshTetrahedronMesher<decimal> aTetrahedronMesher;

    aTetrahedronMesher.generate(aSurfaceMesh, 999, d, d * 0.1, d * 10.0, 1E-3);

    CMshMesh<decimal> aVolumeMesh;
    aVolumeMesh = aTetrahedronMesher.mesh();

    T.setMesh(aVolumeMesh);
    aPosGmsh.save(T, "tetra_volume2.msh", "tetras");

    EXPECT_EQ(95, aVolumeMesh.nbElements());

}

TEST_F(CTestMshTetrahedronMesher, mesh3) {

    const decimal d = 0.125;

    const Integer nu = 5;
    const Integer nv = 3;
    const Integer nw = 3;    

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(nu * d, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(nu * d, nv * d, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, nv * d, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, nw * d);
    CGeoCoordinate<decimal> aVertex6(nu * d, 0.0, nw * d);
    CGeoCoordinate<decimal> aVertex7(nu * d, nv * d, nw * d);
    CGeoCoordinate<decimal> aVertex8(0.0, nv * d, nw * d);

    CGeoHexahedron<decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CMshBasicMesher<decimal> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    CMshMesh<decimal> aSurfaceMesh = aBasicMesher.mesh().extractBoundary(1E-6);

    aSurfaceMesh.generateFaces(1E-5);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tetra_surface3.msh", "tris");

    CMshTetrahedronMesher<decimal> aTetrahedronMesher;

    aTetrahedronMesher.generate(aSurfaceMesh, 999, d * 1.1, d * 1.0, d * 1.2, 1E-3);

    CMshMesh<decimal> aVolumeMesh;
    aVolumeMesh = aTetrahedronMesher.mesh();

    T.setMesh(aVolumeMesh);
    aPosGmsh.save(T, "tetra_volume3.msh", "tetras");

    EXPECT_EQ(281, aVolumeMesh.nbElements());

}

