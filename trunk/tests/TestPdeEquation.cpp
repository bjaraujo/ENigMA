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
#include "PdeEquation.hpp"
#include "PdeField.hpp"

using namespace ENigMA::pde;

class CTestPdeEquation : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestPdeEquation, femSteadyLaplaceLine)
{

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<decimal> aLine(aPoint1, aPoint2);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 10;

    aBasicMesher.generate(aLine, nu);

    EXPECT_EQ(nu, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    // Set BC
    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

    }

    // Steady-state conduction in a line
    CPdeEquation<decimal> aPdeEquation(laplacian<decimal>(T) = 0);

    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i) {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        EXPECT_NEAR(aNode.x(), T.u(i), 1E-2);
    }

}

TEST_F(CTestPdeEquation, fdmSteadyLaplaceLine)
{

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<decimal> aLine(aPoint1, aPoint2);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 10;

    aBasicMesher.generate(aLine, nu);

    EXPECT_EQ(nu, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FDM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    // Set BC
    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

    }

    // Steady-state conduction in a line
    CPdeEquation<decimal> aPdeEquation(laplacian<decimal>(T) = 0);

    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i) {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        EXPECT_NEAR(aNode.x(), T.u(i), 1E-1);
    }

}

TEST_F(CTestPdeEquation, fvmSteadyLaplaceLine)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.1, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.1, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, 0.1);
    CGeoCoordinate<decimal> aVertex6(1.0, 0.0, 0.1);
    CGeoCoordinate<decimal> aVertex7(1.0, 0.1, 0.1);
    CGeoCoordinate<decimal> aVertex8(0.0, 0.1, 0.1);

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

    const Integer nu = 10;
    const Integer nv = 1;
    const Integer nw = 1;

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    aBasicMesher.mesh().generateFaces(1E-12);

    aBasicMesher.mesh().calculateFaceCentroid();
    aBasicMesher.mesh().calculateElementCentroid();

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FVM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_ELEMENT_CENTER);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC and initial conditions
    T.u.resize(T.mesh().nbElements());

    for (Integer i = 0; i < T.mesh().nbElements(); ++i)
    {
        T.u(i) = 0.0;
    }

    for (Integer i = 0; i < T.mesh().nbFaces(); ++i)
    {

        Integer aFaceId = T.mesh().faceId(i);

        if (fabs(T.mesh().faceCentroid(aFaceId).x() - 0.0) < 1E-6)
        {
            CPdeBoundaryCondition<decimal> aFixedTemperature(BT_GENERIC_FIXED_VALUE);
            aFixedTemperature.addCondition(CT_GENERIC_FIXED_VALUE, 0.0);
            T.addBCFace(aFaceId, aFixedTemperature);
        }

        if (fabs(T.mesh().faceCentroid(aFaceId).x() - 1.0) < 1E-6)
        {
            CPdeBoundaryCondition<decimal> aFixedTemperature(BT_GENERIC_FIXED_VALUE);
            aFixedTemperature.addCondition(CT_GENERIC_FIXED_VALUE, 1.0);
            T.addBCFace(aFaceId, aFixedTemperature);
        }

    }

    // Steady-state conduction in a line
    CPdeEquation<decimal> aPdeEquation(laplacian<decimal>(T) = 0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    for (Integer i = 0; i < T.mesh().nbElements(); ++i) {

        Integer anElementId = T.mesh().elementId(i);

        CGeoCoordinate<decimal> aCentroid = T.mesh().elementCentroid(anElementId);

        EXPECT_NEAR(aCentroid.x(), T.u(i), 1E-2);
    }

}

TEST_F(CTestPdeEquation, bemSteadyLaplaceLine)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.1);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.1);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.1, 0.1);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.1, 0.1);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex6(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex7(1.0, 0.1, 0.0);
    CGeoCoordinate<decimal> aVertex8(0.0, 0.1, 0.0);

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

    const Integer nu = 10;
    const Integer nv = 1;
    const Integer nw = 1;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    CMshMesh<decimal> aSurfaceMesh = aBasicMesher.mesh().extractBoundary(1E-12);
    aSurfaceMesh.invert();

    CPdeField<decimal> T;

    T.setMesh(aSurfaceMesh);
    T.setDiscretMethod(DM_BEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_ELEMENT_CENTER);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC
    T.mesh().calculateElementCentroid();

    for (Integer i = 0; i < T.mesh().nbElements(); ++i) {

        Integer anElementId = T.mesh().elementId(i);

        CGeoCoordinate<decimal> aCentroid = T.mesh().elementCentroid(anElementId);

        if (fabs(aCentroid.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aCentroid.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);
    }

    // Steady-state conduction in a box
    CPdeEquation<decimal> aPdeEquation(laplacian<decimal>(T) = 0);

    aPdeEquation.solve(T);

    for (Integer i = 0; i < T.mesh().nbElements(); ++i) {

        Integer anElementId = T.mesh().elementId(i);

        CGeoCoordinate<decimal> aCentroid = T.mesh().elementCentroid(anElementId);

        EXPECT_NEAR(aCentroid.x(), T.u(i), 1E-2);
    }

}

TEST_F(CTestPdeEquation, femTransientLaplaceLine)
{

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<decimal> aLine(aPoint1, aPoint2);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 10;

    aBasicMesher.generate(aLine, nu);

    EXPECT_EQ(nu, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    // Set BC and initial temperature

    T.u.resize(T.mesh().nbNodes());

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

        T.u(i) = 0.0;

    }

    // Transient conduction in steady state
    CPdeEquation<decimal> aPdeEquation(static_cast<decimal>(1 / 10.0) * ddt<decimal>(T) - laplacian<decimal>(T) = 0);

    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i) {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        EXPECT_NEAR(aNode.x(), T.u(i), 1E-2);
    }

}

TEST_F(CTestPdeEquation, fdmTransientLaplaceLine)
{

    CGeoCoordinate<decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<decimal> aLine(aPoint1, aPoint2);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 10;

    aBasicMesher.generate(aLine, nu);

    EXPECT_EQ(nu, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FDM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    // Set BC and initial temperature

    T.u.resize(T.mesh().nbNodes());

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

        T.u(i) = 0.0;

    }

    // Transient conduction in steady state
    CPdeEquation<decimal> aPdeEquation(static_cast<decimal>(1 / 10.0) * ddt<decimal>(T) - laplacian<decimal>(T) = 0);

    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i) {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        EXPECT_NEAR(aNode.x(), T.u(i), 1E-1);
    }

}

TEST_F(CTestPdeEquation, femSteadyLaplaceRectangle1)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.5, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.5, 0.0);

    CGeoQuadrilateral<decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;

    aBasicMesher.generate(aQuadrilateral, nu, nv);

    EXPECT_EQ(nu * nv, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    // Set BC
    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

    }

    // Steady-state conduction in a rectangular plate
    CPdeEquation<decimal> aPdeEquation(laplacian<decimal>(T) = 0);

    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    EXPECT_NEAR(0.00, T.u(4), 1E-2);
    EXPECT_NEAR(0.33, T.u(5), 1E-2);
    EXPECT_NEAR(0.66, T.u(6), 1E-2);
    EXPECT_NEAR(1.00, T.u(7), 1E-2);

}

TEST_F(CTestPdeEquation, femSteadyLaplaceRectangle2)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.5, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.5, 0.0);

    CGeoQuadrilateral<decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;

    aBasicMesher.generate(aQuadrilateral, nu, nv, true);

    EXPECT_EQ(nu * nv * 2, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    // Set BC
    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

    }

    // Steady-state conduction in a rectangular plate
    CPdeEquation<decimal> aPdeEquation(laplacian<decimal>(T) = 0);

    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    EXPECT_NEAR(0.00, T.u(4), 1E-2);
    EXPECT_NEAR(0.33, T.u(5), 1E-2);
    EXPECT_NEAR(0.66, T.u(6), 1E-2);
    EXPECT_NEAR(1.00, T.u(7), 1E-2);

}

TEST_F(CTestPdeEquation, femSteadyLaplaceBox1)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.5, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.5, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, 0.5);
    CGeoCoordinate<decimal> aVertex6(1.0, 0.0, 0.5);
    CGeoCoordinate<decimal> aVertex7(1.0, 0.5, 0.5);
    CGeoCoordinate<decimal> aVertex8(0.0, 0.5, 0.5);

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

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 2;

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    EXPECT_EQ(nu * nv * nw, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    // Set BC
    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

    }

    // Steady-state conduction in a box
    CPdeEquation<decimal> aPdeEquation(laplacian<decimal>(T) = 0);

    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    EXPECT_NEAR(0.00, T.u(4), 1E-2);
    EXPECT_NEAR(0.33, T.u(5), 1E-2);
    EXPECT_NEAR(0.66, T.u(6), 1E-2);
    EXPECT_NEAR(1.00, T.u(7), 1E-2);

}

TEST_F(CTestPdeEquation, femSteadyLaplaceBox2)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.5, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 0.5, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, 0.5);
    CGeoCoordinate<decimal> aVertex6(1.0, 0.0, 0.5);
    CGeoCoordinate<decimal> aVertex7(1.0, 0.5, 0.5);
    CGeoCoordinate<decimal> aVertex8(0.0, 0.5, 0.5);

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

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 2;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    EXPECT_EQ(nu * nv * nw * 6, aBasicMesher.mesh().nbElements());

    CPdeField<decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    // Set BC
    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        CMshNode<decimal> aNode = T.mesh().node(T.mesh().nodeId(i));

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

    }

    // Steady-state conduction in a box
    CPdeEquation<decimal> aPdeEquation(laplacian<decimal>(T) = 0);

    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    EXPECT_NEAR(0.00, T.u(4), 1E-2);
    EXPECT_NEAR(0.33, T.u(5), 1E-2);
    EXPECT_NEAR(0.66, T.u(6), 1E-2);
    EXPECT_NEAR(1.00, T.u(7), 1E-2);

}
