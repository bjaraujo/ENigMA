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
#include "FvmPisoSolver.hpp"

using namespace ENigMA::fvm;

class CTestFvmPiso : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFvmPiso, hydroPressure) {

    CGeoCoordinate<decimal> aVertex1(+0.00, +0.00, +1.0);
    CGeoCoordinate<decimal> aVertex2(+1.00, +0.00, +1.0);
    CGeoCoordinate<decimal> aVertex3(+1.00, +1.00, +1.0);
    CGeoCoordinate<decimal> aVertex4(+0.00, +1.00, +1.0);
    CGeoCoordinate<decimal> aVertex5(+0.00, +0.00, -1.0);
    CGeoCoordinate<decimal> aVertex6(+1.00, +0.00, -1.0);
    CGeoCoordinate<decimal> aVertex7(+1.00, +1.00, -1.0);
    CGeoCoordinate<decimal> aVertex8(+0.00, +1.00, -1.0);

    CGeoHexahedron<decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);
    
    Integer nx = 10;
    Integer ny = 80;

    CMshBasicMesher<decimal> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nx, ny, 1);
    
    aBasicMesher.mesh().generateFaces(1E-12);

    aBasicMesher.mesh().calculateFaceCentroid();
    aBasicMesher.mesh().calculateElementCentroid();

    CFvmMesh<decimal> aFvmMesh(aBasicMesher.mesh());

    CFvmPisoSolver<decimal> aPisoSolver(aFvmMesh);

    decimal g = -9.8;
    aPisoSolver.setGravity(0.0, g, 0.0);

    decimal mu = 0.1; // dynamic viscosity
    decimal rho = 1000.0; // density
 
    aPisoSolver.setMaterialProperties(rho, mu);

    std::vector<Integer> sFaceIds;

    sFaceIds.clear();
    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<decimal> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if ((aFace.centroid().x() - 0.0) < 1E-3 ||
            (aFace.centroid().x() - 1.0) < 1E-3 ||
            (aFace.centroid().y() - 0.0) < 1E-3 ||
            (aFace.centroid().y() - 1.0) < 1E-3)
        {
            sFaceIds.push_back(aFaceId);
        }

    }

    aPisoSolver.setBoundaryVelocity(sFaceIds, BT_WALL_NO_SLIP, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(sFaceIds, BT_WALL_NO_SLIP, 0.0);

    decimal dt = 1E-3;
    Integer nIter = 10;
    
    for (Integer i = 0; i < nIter; ++i)
    {
        aPisoSolver.iterate(dt);
    }

    decimal p = 0.0;
    
    for (Integer i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
    {
        unsigned int aControlVolumeId = aFvmMesh.controlVolumeId(i);
        p = std::max(p, aPisoSolver.p(aControlVolumeId));
    }
     
    EXPECT_NEAR(rho*std::fabs(g), p, 200);

}

