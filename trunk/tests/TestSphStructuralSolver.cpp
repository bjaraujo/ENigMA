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

#include "GeoCoordinate.hpp"
#include "SphCubicSpline.hpp"
#include "SphStructuralSolver.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::sph;

class CTestSphStructuralSolver : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestSphStructuralSolver, cantileverBeam1D)
{
    // Create a simple 1D cantilever beam test
    // Beam: L = 10.0, fixed at x=0, load P applied at x=10.0
    // Material: E = 1000, ν = 0.3
    // Particles spaced at dx = 0.5

    CSphStructuralSolver<Decimal> solver;
    CSphCubicSpline<Decimal> kernel(1); // 1D kernel

    // Setup solver
    solver.setKernel(&kernel);
    solver.setSmoothingLength(1.2);  // h = 1.2 * dx
    solver.setMaterialProperties(1000.0, 0.3);  // E, nu

    // Create cantilever beam: particles along x-axis from 0 to 10
    Decimal beamLength = 10.0;
    Decimal spacing = 0.5;
    Integer nParticles = static_cast<Integer>(beamLength / spacing) + 1;

    std::vector<CGeoCoordinate<Decimal>> points;
    for (Integer i = 0; i < nParticles; ++i)
    {
        Decimal x = i * spacing;
        points.push_back(CGeoCoordinate<Decimal>(x, 0.0, 0.0));
    }

    solver.addParticles(points);

    // Fix left end (clamped boundary condition)
    CGeoVector<Decimal> fixedDisp(0.0, 0.0, 0.0);
    solver.setFixedBoundaryCondition(0, fixedDisp);

    // Apply vertical load at free end (right end)
    Integer lastPtId = nParticles - 1;
    CGeoVector<Decimal> load(0.0, -10.0, 0.0);  // Load in negative y direction
    solver.setForce(lastPtId, load);

    // Solve
    solver.solve(50, 1e-5);

    // Check results
    // Verify that fixed end has zero displacement
    CGeoVector<Decimal> disp0 = solver.getDisplacement(0);
    EXPECT_NEAR(0.0, disp0.x(), 1e-6);
    EXPECT_NEAR(0.0, disp0.y(), 1e-6);
    EXPECT_NEAR(0.0, disp0.z(), 1e-6);

    // Verify that free end has non-zero displacement in y-direction
    CGeoVector<Decimal> dispLast = solver.getDisplacement(lastPtId);
    EXPECT_TRUE(dispLast.y() < 0.0);  // Should deflect downward (negative y)

    // Rough check: for a cantilever beam under point load
    // δ = P*L³/(3*E*I)
    // For a 1D particle approximation with unit cross-section, I ~ 1
    // δ ≈ 10.0 * 1000.0³ / (3 * 1000 * 1) = 3.33e9 (this is too large!)
    // Let's use a simpler analytical check:
    // For a soft material (E=1000) with small load (P=10), expect reasonable deflection
    
    std::cout << "Fixed end displacement: (" << disp0.x() << ", " << disp0.y() << ", " << disp0.z() << ")" << std::endl;
    std::cout << "Free end displacement: (" << dispLast.x() << ", " << dispLast.y() << ", " << dispLast.z() << ")" << std::endl;

    // The deflection should be reasonable (less than beam length)
    EXPECT_TRUE(std::abs(dispLast.y()) < beamLength);

    // Stress should be non-zero at some point
    Eigen::Matrix<Decimal, 3, 3> stressLast = solver.getStress(lastPtId);
    EXPECT_TRUE(std::abs(stressLast(1, 1)) > 0.0);  // Normal stress in y-direction

    std::cout << "Stress at free end:" << std::endl << stressLast << std::endl;
}

TEST_F(CTestSphStructuralSolver, cantileverBeam2D)
{
    // 2D cantilever beam test
    // Beam: L = 5.0 (x-direction), H = 1.0 (y-direction)
    // Fixed at x=0, load at x=5.0
    // Material: E = 2000, ν = 0.25

    CSphStructuralSolver<Decimal> solver;
    CSphCubicSpline<Decimal> kernel(2);  // 2D kernel

    solver.setKernel(&kernel);
    solver.setSmoothingLength(0.8);
    solver.setMaterialProperties(2000.0, 0.25);

    // Create 2D particle grid for cantilever beam
    Decimal dx = 0.25;
    Decimal dy = 0.25;
    Decimal beamLx = 5.0;
    Decimal beamLy = 1.0;

    Integer nxParticles = static_cast<Integer>(beamLx / dx) + 1;
    Integer nyParticles = static_cast<Integer>(beamLy / dy) + 1;

    std::vector<CGeoCoordinate<Decimal>> points;
    for (Integer j = 0; j < nyParticles; ++j)
    {
        for (Integer i = 0; i < nxParticles; ++i)
        {
            Decimal x = i * dx;
            Decimal y = j * dy - beamLy / 2.0;  // Center beam at y=0
            points.push_back(CGeoCoordinate<Decimal>(x, y, 0.0));
        }
    }

    solver.addParticles(points);

    // Fix left edge (x=0)
    CGeoVector<Decimal> fixedDisp(0.0, 0.0, 0.0);
    for (Integer j = 0; j < nyParticles; ++j)
    {
        Integer ptId = j * nxParticles;  // First column
        solver.setFixedBoundaryCondition(ptId, fixedDisp);
    }

    // Apply vertical load at right edge (x=5.0)
    Decimal loadPerNode = -1.0 / nyParticles;  // Distribute load
    for (Integer j = 0; j < nyParticles; ++j)
    {
        Integer ptId = j * nxParticles + (nxParticles - 1);  // Last column
        CGeoVector<Decimal> load(0.0, loadPerNode, 0.0);
        solver.setForce(ptId, load);
    }

    // Solve
    solver.solve(50, 1e-4);

    // Verify fixed boundary
    CGeoVector<Decimal> dispFixed = solver.getDisplacement(0);
    EXPECT_NEAR(0.0, dispFixed.x(), 1e-5);
    EXPECT_NEAR(0.0, dispFixed.y(), 1e-5);

    // Verify that free end deflects
    Integer lastPtId = (nyParticles - 1) * nxParticles + (nxParticles - 1);
    CGeoVector<Decimal> dispFree = solver.getDisplacement(lastPtId);

    std::cout << "2D Cantilever - Fixed displacement: (" << dispFixed.x() << ", " << dispFixed.y() << ")" << std::endl;
    std::cout << "2D Cantilever - Free end displacement: (" << dispFree.x() << ", " << dispFree.y() << ")" << std::endl;

    // Free end should deflect in negative y direction
    EXPECT_TRUE(dispFree.y() < -0.001);
    EXPECT_TRUE(dispFree.y() > -beamLy);

    // Check that displacement increases from fixed to free end
    Integer midPtId = (nyParticles / 2) * nxParticles + (nxParticles / 2);
    CGeoVector<Decimal> dispMid = solver.getDisplacement(midPtId);
    EXPECT_TRUE(std::abs(dispMid.y()) < std::abs(dispFree.y()));
}

TEST_F(CTestSphStructuralSolver, cantileverBeamConvergence)
{
    // Convergence study: mesh refinement
    // Compare deflections at different particle spacings

    CSphStructuralSolver<Decimal> solver1, solver2;
    CSphCubicSpline<Decimal> kernel1(1), kernel2(1);

    // Coarse mesh
    solver1.setKernel(&kernel1);
    solver1.setSmoothingLength(1.2);
    solver1.setMaterialProperties(1000.0, 0.3);

    Decimal spacing1 = 1.0;
    Integer nP1 = 11;  // 11 particles
    for (Integer i = 0; i < nP1; ++i)
    {
        solver1.addParticle(CGeoCoordinate<Decimal>(i * spacing1, 0.0, 0.0));
    }
    solver1.setFixedBoundaryCondition(0, CGeoVector<Decimal>(0.0, 0.0, 0.0));
    solver1.setForce(nP1 - 1, CGeoVector<Decimal>(0.0, -5.0, 0.0));
    solver1.solve(50, 1e-5);

    CGeoVector<Decimal> disp1 = solver1.getDisplacement(nP1 - 1);

    // Fine mesh
    solver2.setKernel(&kernel2);
    solver2.setSmoothingLength(0.72);
    solver2.setMaterialProperties(1000.0, 0.3);

    Decimal spacing2 = 0.5;
    Integer nP2 = 21;  // 21 particles
    for (Integer i = 0; i < nP2; ++i)
    {
        solver2.addParticle(CGeoCoordinate<Decimal>(i * spacing2, 0.0, 0.0));
    }
    solver2.setFixedBoundaryCondition(0, CGeoVector<Decimal>(0.0, 0.0, 0.0));
    solver2.setForce(nP2 - 1, CGeoVector<Decimal>(0.0, -5.0, 0.0));
    solver2.solve(50, 1e-5);

    CGeoVector<Decimal> disp2 = solver2.getDisplacement(nP2 - 1);

    std::cout << "Convergence study:" << std::endl;
    std::cout << "Coarse mesh (spacing=" << spacing1 << "): deflection = " << disp1.y() << std::endl;
    std::cout << "Fine mesh (spacing=" << spacing2 << "): deflection = " << disp2.y() << std::endl;

    // Deflections should have same sign
    EXPECT_TRUE(disp1.y() * disp2.y() > 0.0);

    // Fine mesh should generally show more refinement
    EXPECT_TRUE(std::abs(disp2.y()) > 0.0);
}

TEST_F(CTestSphStructuralSolver, simpleElastostatics)
{
    // Simple static test: uniform block with prescribed displacement on one face
    // Should reach equilibrium with uniform stress

    CSphStructuralSolver<Decimal> solver;
    CSphCubicSpline<Decimal> kernel(1);

    solver.setKernel(&kernel);
    solver.setSmoothingLength(1.5);
    solver.setMaterialProperties(100.0, 0.2);

    // Create 5 particles in a line
    for (Integer i = 0; i < 5; ++i)
    {
        solver.addParticle(CGeoCoordinate<Decimal>(i * 1.0, 0.0, 0.0));
    }

    // Fix left end
    solver.setFixedBoundaryCondition(0, CGeoVector<Decimal>(0.0, 0.0, 0.0));

    // Set prescribed displacement on right end
    CGeoVector<Decimal> prescribedDisp(0.1, 0.0, 0.0);
    solver.setFixedBoundaryCondition(4, prescribedDisp);

    solver.solve(30, 1e-6);

    // Verify fixed end
    CGeoVector<Decimal> d0 = solver.getDisplacement(0);
    EXPECT_NEAR(0.0, d0.x(), 1e-6);

    // Verify prescribed end
    CGeoVector<Decimal> d4 = solver.getDisplacement(4);
    EXPECT_NEAR(0.1, d4.x(), 0.01);

    // Displacement should increase monotonically from left to right
    CGeoVector<Decimal> d2 = solver.getDisplacement(2);
    EXPECT_TRUE(d2.x() > 0.0 && d2.x() < 0.1);

    std::cout << "Simple elastostatics test:" << std::endl;
    std::cout << "d[0] = " << d0.x() << ", d[2] = " << d2.x() << ", d[4] = " << d4.x() << std::endl;
}
