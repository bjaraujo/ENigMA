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
#include "FemCbsSolver.hpp"

using namespace ENigMA::fem;
using namespace ENigMA::mesh;

class CTestFemCbs : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFemCbs, hydroPressure) {

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

    aBasicMesher.generate(aQuadrilateral, 60, 60, true);
    CMshMesh<decimal> aMesh = aBasicMesher.mesh();
    
    CFemCbsSolver<decimal, 2> aCbsSolver(aMesh);

    decimal g = -9.8;
    aCbsSolver.setGravity(0.0, -9.8);

    decimal mu = 0.1; // dynamic viscosity
    decimal rho = 1000.0; // density
 
    aCbsSolver.setMaterialProperties(rho, mu);

    
     
     
    //EXPECT_NEAR(rho*fabs(g), p, 200);

}

