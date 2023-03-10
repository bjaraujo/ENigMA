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

using namespace ENigMA::mesh;

class CTestMshBasicMesher : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestMshBasicMesher, meshLine) {

    CGeoCoordinate<Decimal> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<Decimal> aLine(aPoint1, aPoint2);

    aLine.calculateLength();

    EXPECT_NEAR(1.0, aLine.length(), 1E-12);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer ne = 20;

    aBasicMesher.generate(aLine, ne);

    EXPECT_EQ(ne, aBasicMesher.mesh().nbElements());

    Decimal sum_length = 0.0;

    for (Integer i = 0; i < aBasicMesher.mesh().nbElements(); ++i)
    {

        Integer elementId = aBasicMesher.mesh().elementId(i);

        CMshElement<Decimal> aElement = aBasicMesher.mesh().element(elementId);

        CGeoLine<Decimal> aLine;

        for (Integer j = 0; j < aElement.nbNodeIds(); ++j)
        {

            CMshNode<Decimal> aNode = aBasicMesher.mesh().node(aElement.nodeId(j));
            CGeoCoordinate<Decimal> aPoint = aNode;
            if (j == 0)
                aLine.setStartPoint(aPoint);
            else
                aLine.setEndPoint(aPoint);
        }

        aLine.calculateLength();

        sum_length += aLine.length();

    }

    EXPECT_NEAR(1.0, sum_length, 1E-12);

}

TEST_F(CTestMshBasicMesher, meshQuadrilateral1) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoQuadrilateral<Decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    aQuadrilateral.calculateArea();

    EXPECT_NEAR(1.0, aQuadrilateral.area(), 1E-12);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;

    aBasicMesher.generate(aQuadrilateral, nu, nv);

    EXPECT_EQ(nu * nv, aBasicMesher.mesh().nbElements());

    Decimal sum_area = 0.0;

    for (Integer i = 0; i < aBasicMesher.mesh().nbElements(); ++i)
    {

        Integer elementId = aBasicMesher.mesh().elementId(i);

        CMshElement<Decimal> aElement = aBasicMesher.mesh().element(elementId);

        CGeoQuadrilateral<Decimal> aQuadrilateral;

        for (Integer j = 0; j < aElement.nbNodeIds(); ++j)
        {

            CMshNode<Decimal> aNode = aBasicMesher.mesh().node(aElement.nodeId(j));
            CGeoCoordinate<Decimal> aVertex = aNode;
            aQuadrilateral.addVertex(aVertex);

        }

        aQuadrilateral.calculateArea();

        sum_area += aQuadrilateral.area();

    }

    EXPECT_NEAR(1.0, sum_area, 1E-12);

}

TEST_F(CTestMshBasicMesher, meshQuadrilateral2) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);

    CGeoQuadrilateral<Decimal> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    aQuadrilateral.calculateArea();

    EXPECT_NEAR(1.0, aQuadrilateral.area(), 1E-6);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;

    aBasicMesher.generate(aQuadrilateral, nu, nv, true);

    EXPECT_EQ(nu * nv * 2, aBasicMesher.mesh().nbElements());

    Decimal sum_area = 0.0;

    for (Integer i = 0; i < aBasicMesher.mesh().nbElements(); ++i)
    {

        Integer elementId = aBasicMesher.mesh().elementId(i);

        CMshElement<Decimal> aElement = aBasicMesher.mesh().element(elementId);

        CGeoTriangle<Decimal> aTriangle;

        for (Integer j = 0; j < aElement.nbNodeIds(); ++j)
        {

            CMshNode<Decimal> aNode = aBasicMesher.mesh().node(aElement.nodeId(j));
            CGeoCoordinate<Decimal> aVertex = aNode;
            aTriangle.addVertex(aVertex);

        }

        aTriangle.calculateArea();

        sum_area += aTriangle.area();

    }

    EXPECT_NEAR(1.0, sum_area, 1E-6);

}

TEST_F(CTestMshBasicMesher, meshHexahedron1) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex7(1.0, 1.0, -1.0);
    CGeoCoordinate<Decimal> aVertex8(0.0, 1.0, -1.0);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    aHexahedron.calculateVolume();

    EXPECT_NEAR(1.0, aHexahedron.volume(), 1E-6);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 2;

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    EXPECT_EQ(nu * nv * nw, aBasicMesher.mesh().nbElements());

    aBasicMesher.mesh().generateFaces(1E-6);

    EXPECT_EQ(2 * (nu * nv + nv * nw + nw * nu), aBasicMesher.mesh().nbBoundaryFaces());

    Decimal sum_volume = 0.0;

    for (Integer i = 0; i < aBasicMesher.mesh().nbElements(); ++i)
    {

        Integer elementId = aBasicMesher.mesh().elementId(i);

        CMshElement<Decimal> aElement = aBasicMesher.mesh().element(elementId);

        CGeoHexahedron<Decimal> aHexahedron;

        for (Integer j = 0; j < aElement.nbNodeIds(); ++j)
        {

            CMshNode<Decimal> aNode = aBasicMesher.mesh().node(aElement.nodeId(j));
            CGeoCoordinate<Decimal> aVertex = aNode;
            aHexahedron.addVertex(aVertex);

        }

        aHexahedron.calculateVolume();

        sum_volume += aHexahedron.volume();

    }

    EXPECT_NEAR(1.0, sum_volume, 1E-6);

}

TEST_F(CTestMshBasicMesher, meshHexahedron2) {

    CGeoCoordinate<Decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<Decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<Decimal> aVertex5(0.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex6(1.0, 0.0, -1.0);
    CGeoCoordinate<Decimal> aVertex7(1.0, 1.0, -1.0);
    CGeoCoordinate<Decimal> aVertex8(0.0, 1.0, -1.0);

    CGeoHexahedron<Decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    aHexahedron.calculateVolume();

    EXPECT_NEAR(1.0, aHexahedron.volume(), 1E-6);

    CMshBasicMesher<Decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 2;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    EXPECT_EQ(nu * nv * nw * 6, aBasicMesher.mesh().nbElements());

    Decimal sum_volume = 0.0;

    for (Integer i = 0; i < aBasicMesher.mesh().nbElements(); ++i)
    {

        Integer elementId = aBasicMesher.mesh().elementId(i);

        CMshElement<Decimal> aElement = aBasicMesher.mesh().element(elementId);

        CGeoTetrahedron<Decimal> aTetrahedron;

        for (Integer j = 0; j < aElement.nbNodeIds(); ++j)
        {

            CMshNode<Decimal> aNode = aBasicMesher.mesh().node(aElement.nodeId(j));
            CGeoCoordinate<Decimal> aVertex = aNode;
            aTetrahedron.addVertex(aVertex);

        }

        aTetrahedron.calculateVolume();

        sum_volume += aTetrahedron.volume();

    }

    EXPECT_NEAR(1.0, sum_volume, 1E-6);

}
