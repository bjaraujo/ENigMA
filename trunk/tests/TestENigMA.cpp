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

class CTestENigMA : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestENigMA, update) {

    /*
    CGeoPointd aPoint1(0.0, 0.0, 0.0);
    CGeoPointd aPoint2(1.0, 0.0, 0.0);

    CGeoLined aLine(&aPoint1, &aPoint2);

    CMshBasicMesherd aBasicMesher;

    const Integer ne = 3;

    CMshMeshd *aMesh = aBasicMesher.meshGeometry(aLine, ne);

    EXPECT_EQ(ne, aMesh->nbElements());

    // mesh, 1 Dimension per unknown or DOF, discretisation method
    CPdeField T(aMesh, 1, pde::FEM);

    CPdeEquationd aPdeEquation;
    aPdeEquation.setEquation(pde::ddt(T) + pde::laplacian(T), 0)
    aPdeEquation.setSolver(GMRES);
    aPdeEquation.setSolverTolerance(1E-6);
    aPdeEquation.setSolverIterations(100);

    CDomaind aDomain;

    //aDomain.setVolumeMesh(aMesh);
    //aDomain.setBoundaryMesh(aBndMesh);
    //aDomain.extractBoundary();
    aDomain.addTimeInterval(dt);
    aDomain.addMaterial(material1);
    aDomain.addMaterial(material2);
    aDomain.addPDE(pde1);
    aDomain.addPDE(pde2);
    aDomain.solve();

    CSolution aSolution
    aSolution.addDomain(aDomain);
    aSolution.solve();
    */

}

