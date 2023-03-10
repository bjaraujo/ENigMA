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
#include "PdeField.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::pde;

class CTestPdeGeometricField : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestPdeGeometricField, set) {

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<Decimal> aLine(aPoint1, aPoint2);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer ne = 3;

    aBasicMesher.generate(aLine, ne);

    EXPECT_EQ(ne, aBasicMesher.mesh().nbElements());

    CPdeField<Decimal> T;

    T.setMesh(aBasicMesher.mesh());
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    EXPECT_EQ(ne, T.mesh().nbElements());

}

