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

#include "PosGmsh.hpp"

using namespace ENigMA::fem;
using namespace ENigMA::mesh;

using namespace ENigMA::post;

class CTestFemCbs : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFemCbs, hydroPressure) {

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

    aBasicMesher.generate(aQuadrilateral, 10, 40, true);
    
    CFemCbsSolver<Decimal, 2> aCbsSolver(aBasicMesher.mesh());

    Decimal g = -9.8;
    aCbsSolver.setGravity(0.0, -9.8);

    Decimal mu = 0.1; // dynamic viscosity
    Decimal rho = 1000.0; // density
 
    aCbsSolver.setMaterialProperties(rho, mu);

    for (Integer i = 0; i < aBasicMesher.mesh().nbNodes(); ++i)
    {
        Integer aNodeId = aBasicMesher.mesh().nodeId(i);
        CMshNode<Decimal> aNode = aBasicMesher.mesh().node(aNodeId);

        if (std::fabs(aNode.x() - 0.0) < 1E-3 ||
            std::fabs(aNode.x() - 1.0) < 1E-3)
        {
            aCbsSolver.u().setFixedValue(i, 0.0);
        }

        if (std::fabs(aNode.y() - 0.0) < 1E-3)
        {
            aCbsSolver.u().setFixedValue(i, 0.0);
            aCbsSolver.v().setFixedValue(i, 0.0);
        }

        if (std::fabs(aNode.y() - 1.0) < 1E-3)
        {
            aCbsSolver.u().setFixedValue(i, 0.0);
            aCbsSolver.v().setFixedValue(i, 0.0);
            aCbsSolver.p().setFixedValue(i, 0.0);
        }

        aCbsSolver.u().setValue(i, 0.0);
        aCbsSolver.v().setValue(i, 0.0);

        aCbsSolver.p().setValue(i, 0.0);
    }

    Decimal dt = 1E-3;
    Integer nIter = 2;

    for (Integer i = 0; i < nIter; ++i)
    {
        aCbsSolver.iterate(dt);
    }

    Decimal p = 0.0;

    for (Integer i = 0; i < aBasicMesher.mesh().nbNodes(); ++i)
    {
        p = std::max(p, aCbsSolver.p().value(i));
    }

    EXPECT_NEAR(rho*std::fabs(g), p, 1000);
}

