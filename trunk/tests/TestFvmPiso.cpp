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
#include "FvmTemperatureSolver.hpp"

using namespace ENigMA::fvm;

class CTestFvmPiso : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFvmPiso, hydroPressure) {

    CGeoCoordinate<Decimal> aVertex1(+0.00, +0.00, +1.0);
    CGeoCoordinate<Decimal> aVertex2(+1.00, +0.00, +1.0);
    CGeoCoordinate<Decimal> aVertex3(+1.00, +1.00, +1.0);
    CGeoCoordinate<Decimal> aVertex4(+0.00, +1.00, +1.0);
    CGeoCoordinate<Decimal> aVertex5(+0.00, +0.00, -1.0);
    CGeoCoordinate<Decimal> aVertex6(+1.00, +0.00, -1.0);
    CGeoCoordinate<Decimal> aVertex7(+1.00, +1.00, -1.0);
    CGeoCoordinate<Decimal> aVertex8(+0.00, +1.00, -1.0);

    CGeoHexahedron<Decimal> aHexahedron;

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

    CMshBasicMesher<Decimal> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nx, ny, 1);
    
    aBasicMesher.mesh().generateFaces(1E-12);

    aBasicMesher.mesh().calculateFaceCentroid();
    aBasicMesher.mesh().calculateElementCentroid();

    CFvmMesh<Decimal> aFvmMesh(aBasicMesher.mesh());

    CFvmPisoSolver<Decimal> aPisoSolver(aFvmMesh);

    Decimal g = -9.8;
    aPisoSolver.setGravity(0.0, g, 0.0);

    Decimal mu = 0.1; // dynamic viscosity
    Decimal rho = 1000.0; // density
 
    aPisoSolver.setMaterialProperties(rho, mu);

    std::vector<Integer> sFaceIds;

    sFaceIds.clear();
    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<Decimal> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if ((aFace.centroid().x() - 0.0) < 1E-3 ||
            (aFace.centroid().x() - 1.0) < 1E-3 ||
            (aFace.centroid().y() - 0.0) < 1E-3 ||
            (aFace.centroid().y() - 1.0) < 1E-3)
        {
            sFaceIds.emplace_back(aFaceId);
        }

    }

    aPisoSolver.setBoundaryVelocity(sFaceIds, BT_WALL_NO_SLIP, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(sFaceIds, BT_WALL_NO_SLIP, 0.0);

    Decimal dt = 1E-3;
    Integer nIter = 10;
    
    for (Integer i = 0; i < nIter; ++i)
    {
        aPisoSolver.iterate(dt);
    }

    Decimal p = 0.0;
    
    for (Integer i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
    {
        unsigned int aControlVolumeId = aFvmMesh.controlVolumeId(i);
        p = std::max(p, aPisoSolver.p(aControlVolumeId));
    }
     
    EXPECT_NEAR(rho*std::fabs(g), p, 200);

}

TEST_F(CTestFvmPiso, channelPressure) {

    Decimal L = 0.1;
    Decimal a = 0.01;

    CGeoCoordinate<Decimal> aVertex1(+0.00, +0.00, +L/2);
    CGeoCoordinate<Decimal> aVertex2(+a, +0.00, +L/2);
    CGeoCoordinate<Decimal> aVertex3(+a, +a, +L/2);
    CGeoCoordinate<Decimal> aVertex4(+0.00, +a, +L/2);
    CGeoCoordinate<Decimal> aVertex5(+0.00, +0.00, -L/2);
    CGeoCoordinate<Decimal> aVertex6(+a, +0.00, -L/2);
    CGeoCoordinate<Decimal> aVertex7(+a, +a, -L/2);
    CGeoCoordinate<Decimal> aVertex8(+0.00, +a, -L/2);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);
    
    Integer nx = 4;
    Integer ny = 4;
    Integer nz = 40;

    CMshBasicMesher<Decimal> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nx, ny, nz);
    
    aBasicMesher.mesh().generateFaces(1E-3);

    aBasicMesher.mesh().calculateFaceCentroid();
    aBasicMesher.mesh().calculateElementCentroid();

    CFvmMesh<Decimal> aFvmMesh(aBasicMesher.mesh());

    CFvmPisoSolver<Decimal> aPisoSolver(aFvmMesh);

    Decimal mu = 1E-4; // dynamic viscosity
    Decimal rho = 1.0; // density
 
    aPisoSolver.setGravity(0.0, 0.0, 0.0);
    aPisoSolver.setMaterialProperties(rho, mu);

    std::vector<Integer> aInletFaceIds;
    std::vector<Integer> aOutletFaceIds;
    std::vector<Integer> aWallFaceIds;

    aInletFaceIds.clear();
    aOutletFaceIds.clear();
    aWallFaceIds.clear();

    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {
        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<Decimal> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (std::abs(aFace.centroid().z() + L/2) < 1E-3)
        {
            aInletFaceIds.emplace_back(aFaceId);
        }
        else if (std::abs(aFace.centroid().z() - L/2) < 1E-3)
        {
            aOutletFaceIds.emplace_back(aFaceId);
        }
        else
        {
            aWallFaceIds.emplace_back(aFaceId);
        }
    }

    Decimal v = 0.01;
    aPisoSolver.setBoundaryVelocity(aInletFaceIds, BT_INLET_FLOW, 0.0, 0.0, v);
    aPisoSolver.setBoundaryPressure(aInletFaceIds, BT_INLET_FLOW, 0.0);
    aPisoSolver.setBoundaryVelocity(aOutletFaceIds, BT_OUTLET, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(aOutletFaceIds, BT_OUTLET, 0.0);
    aPisoSolver.setBoundaryVelocity(aWallFaceIds, BT_WALL_NO_SLIP, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(aWallFaceIds, BT_WALL_NO_SLIP, 0.0);

    Decimal dt = 1E-2;
    Integer nIter = 20;
    
    for (Integer i = 0; i < nIter; ++i)
    {
        aPisoSolver.iterate(dt);
    }

    Decimal pMax = 0.0;
    
    for (Integer i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
    {
        unsigned int aControlVolumeId = aFvmMesh.controlVolumeId(i);
        pMax = std::max(pMax, aPisoSolver.p(aControlVolumeId));
    }

    Decimal pExpected = 28.45 * mu * L * v / (a * a);
    EXPECT_NEAR(pExpected, pMax, 0.15);

}

TEST_F(CTestFvmPiso, channelTemperature) {

    GTEST_SKIP();
    
    Decimal L = 0.1;
    Decimal a = 0.01;

    CGeoCoordinate<Decimal> aVertex1(+0.00, +0.00, +L/2);
    CGeoCoordinate<Decimal> aVertex2(+a, +0.00, +L/2);
    CGeoCoordinate<Decimal> aVertex3(+a, +a, +L/2);
    CGeoCoordinate<Decimal> aVertex4(+0.00, +a, +L/2);
    CGeoCoordinate<Decimal> aVertex5(+0.00, +0.00, -L/2);
    CGeoCoordinate<Decimal> aVertex6(+a, +0.00, -L/2);
    CGeoCoordinate<Decimal> aVertex7(+a, +a, -L/2);
    CGeoCoordinate<Decimal> aVertex8(+0.00, +a, -L/2);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);
    
    Integer nx = 4;
    Integer ny = 4;
    Integer nz = 41;

    CMshBasicMesher<Decimal> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nx, ny, nz);
    
    aBasicMesher.mesh().generateFaces(1E-3);

    aBasicMesher.mesh().calculateFaceCentroid();
    aBasicMesher.mesh().calculateElementCentroid();

    CFvmMesh<Decimal> aFvmMesh(aBasicMesher.mesh());

    CFvmTemperatureSolver<Decimal> aPisoSolver(aFvmMesh);

    Decimal mu = 1E-4; // dynamic viscosity
    Decimal rho = 1.0; // density
    Decimal thcond = 1.0; // thermal conductivity
    Decimal spheat = 1.0; // specific heat
 
    aPisoSolver.setGravity(0.0, 0.0, 0.0);
    aPisoSolver.setMaterialProperties(rho, mu, thcond, spheat);

    std::vector<Integer> aInletFaceIds;
    std::vector<Integer> aOutletFaceIds;
    std::vector<Integer> aWallFaceIds;

    aInletFaceIds.clear();
    aOutletFaceIds.clear();
    aWallFaceIds.clear();

    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<Decimal> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (std::abs(aFace.centroid().z() + L/2) < 1E-3)
        {
            aInletFaceIds.emplace_back(aFaceId);
        }
        else if (std::abs(aFace.centroid().z() - L/2) < 1E-3)
        {
            aOutletFaceIds.emplace_back(aFaceId);
        }
        else
        {
            aWallFaceIds.emplace_back(aFaceId);
        }

    }

    Decimal v = 0.01;
    aPisoSolver.setBoundaryVelocity(aInletFaceIds, BT_INLET_FLOW, 0.0, 0.0, v);
    aPisoSolver.setBoundaryPressure(aInletFaceIds, BT_INLET_FLOW, 0.0);
    aPisoSolver.setBoundaryTemperature(aInletFaceIds, BT_INLET_FLOW, 1.0);
    aPisoSolver.setBoundaryVelocity(aOutletFaceIds, BT_OUTLET, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(aOutletFaceIds, BT_OUTLET, 0.0);
    aPisoSolver.setBoundaryTemperature(aOutletFaceIds, BT_OUTLET, 0.0);
    aPisoSolver.setBoundaryVelocity(aWallFaceIds, BT_WALL_NO_SLIP, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(aWallFaceIds, BT_WALL_NO_SLIP, 0.0);

    Decimal dt = 1E-2;
    Integer nIter = 20;
    
    for (Integer i = 0; i < nIter; ++i)
    {
        aPisoSolver.iterate(dt);
    }

    Decimal TMid = 0.0;
    
    for (Integer i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
    {
        unsigned int aControlVolumeId = aFvmMesh.controlVolumeId(i);
        CGeoCoordinate<Decimal> aCentroid = aFvmMesh.controlVolume(aControlVolumeId).centroid();
        
        if (std::abs(aCentroid.z()) < 1E-3)
        {
            TMid = aPisoSolver.T(aControlVolumeId);
            break;
        }
    }
    
    EXPECT_NEAR(0.5, TMid, 0.1);

}

