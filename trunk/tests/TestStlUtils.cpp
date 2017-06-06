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

#include "TypeDef.hpp"

#include "MshBasicMesher.hpp"
#include "StlUtils.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::stl;
using namespace ENigMA::post;
using namespace ENigMA::pde;

class CTestStlUtils : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestStlUtils, io)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex6(1.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex7(1.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aVertex8(0.0, 1.0, 1.0);

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

    const Integer nu = 2;
    const Integer nv = 2;
    const Integer nw = 2;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    CMshMesh<decimal> aSurfaceMesh = aBasicMesher.mesh().extractBoundary(1E-12);

    CStlUtils<decimal> aStl(aSurfaceMesh);

    aStl.stlFile().stats.type = FT_BINARY;

    EXPECT_TRUE(aStl.save("box.stl"));
    EXPECT_TRUE(aStl.load("box.stl"));

}

TEST_F(CTestStlUtils, connectivity)
{

    CStlUtils<decimal> aStl;

    CGeoCoordinate<decimal> aVertex0(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex1(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(2.0, 1.0, 0.0);

    CStlFacet<decimal> aFacet0;

    aFacet0.addVertex(aVertex0);
    aFacet0.addVertex(aVertex1);
    aFacet0.addVertex(aVertex2);

    aStl.addFacet(0, aFacet0);

    CStlFacet<decimal> aFacet1;

    aFacet1.addVertex(aVertex0);
    aFacet1.addVertex(aVertex2);
    aFacet1.addVertex(aVertex3);

    aStl.addFacet(1, aFacet1);

    CStlFacet<decimal> aFacet2;

    aFacet2.addVertex(aVertex1);
    aFacet2.addVertex(aVertex2);
    aFacet2.addVertex(aVertex5);

    aStl.addFacet(2, aFacet2);

    CStlFacet<decimal> aFacet3;

    aFacet3.addVertex(aVertex1);
    aFacet3.addVertex(aVertex4);
    aFacet3.addVertex(aVertex5);

    aStl.addFacet(3, aFacet3);

    aStl.generateConnectivity(1E-12);

    EXPECT_EQ(-1, aStl.stlFile().facet(0).edge(0).neighbor());
    EXPECT_EQ(2, aStl.stlFile().facet(0).edge(1).neighbor());
    EXPECT_EQ(1, aStl.stlFile().facet(0).edge(2).neighbor());

    EXPECT_EQ(-1, aStl.stlFile().facet(0).edge(0).whichVertexNot());
    EXPECT_EQ(2, aStl.stlFile().facet(0).edge(1).whichVertexNot());
    EXPECT_EQ(2, aStl.stlFile().facet(0).edge(2).whichVertexNot());

    EXPECT_EQ(0, aStl.stlFile().facet(1).edge(0).neighbor());
    EXPECT_EQ(-1, aStl.stlFile().facet(1).edge(1).neighbor());
    EXPECT_EQ(-1, aStl.stlFile().facet(1).edge(2).neighbor());

    EXPECT_EQ(1, aStl.stlFile().facet(1).edge(0).whichVertexNot());
    EXPECT_EQ(-1, aStl.stlFile().facet(1).edge(1).whichVertexNot());
    EXPECT_EQ(-1, aStl.stlFile().facet(1).edge(2).whichVertexNot());

    EXPECT_EQ(0, aStl.stlFile().facet(2).edge(0).neighbor());
    EXPECT_EQ(-1, aStl.stlFile().facet(2).edge(1).neighbor());
    EXPECT_EQ(3, aStl.stlFile().facet(2).edge(2).neighbor());

    EXPECT_EQ(0, aStl.stlFile().facet(2).edge(0).whichVertexNot());
    EXPECT_EQ(-1, aStl.stlFile().facet(2).edge(1).whichVertexNot());
    EXPECT_EQ(1, aStl.stlFile().facet(2).edge(2).whichVertexNot());

    EXPECT_EQ(-1, aStl.stlFile().facet(3).edge(0).neighbor());
    EXPECT_EQ(-1, aStl.stlFile().facet(3).edge(1).neighbor());
    EXPECT_EQ(2, aStl.stlFile().facet(3).edge(2).neighbor());

    EXPECT_EQ(-1, aStl.stlFile().facet(3).edge(0).whichVertexNot());
    EXPECT_EQ(-1, aStl.stlFile().facet(3).edge(1).whichVertexNot());
    EXPECT_EQ(1, aStl.stlFile().facet(3).edge(2).whichVertexNot());

    EXPECT_TRUE(aStl.save("plane_con.stl"));

}

TEST_F(CTestStlUtils, mesh)
{

    CStlUtils<decimal> aStl;

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex6(2.0, 1.0, 0.0);

    CStlFacet<decimal> aFacet1;

    aFacet1.addVertex(aVertex1);
    aFacet1.addVertex(aVertex2);
    aFacet1.addVertex(aVertex3);

    aStl.addFacet(0, aFacet1);

    CStlFacet<decimal> aFacet2;

    aFacet2.addVertex(aVertex1);
    aFacet2.addVertex(aVertex3);
    aFacet2.addVertex(aVertex4);

    aStl.addFacet(1, aFacet2);

    CStlFacet<decimal> aFacet3;

    aFacet3.addVertex(aVertex2);
    aFacet3.addVertex(aVertex6);
    aFacet3.addVertex(aVertex3);

    aStl.addFacet(2, aFacet3);

    CStlFacet<decimal> aFacet4;

    aFacet4.addVertex(aVertex2);
    aFacet4.addVertex(aVertex5);
    aFacet4.addVertex(aVertex6);

    aStl.addFacet(3, aFacet4);

    aStl.generateConnectivity(1E-12);

    CMshMesh<decimal> aSurfaceMesh = aStl.mesh();

    CStlUtils<decimal> aNewStlFile(aSurfaceMesh);

    EXPECT_TRUE(aNewStlFile.save("plane.stl"));

    EXPECT_EQ(6, aSurfaceMesh.nbNodes());
    EXPECT_EQ(4, aSurfaceMesh.nbElements());

}

TEST_F(CTestStlUtils, split)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex6(1.0, 0.0, 1.0);
    CGeoCoordinate<decimal> aVertex7(1.0, 1.0, 1.0);
    CGeoCoordinate<decimal> aVertex8(0.0, 1.0, 1.0);

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

    const Integer nu = 2;
    const Integer nv = 2;
    const Integer nw = 2;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    CMshMesh<decimal> aSurfaceMesh = aBasicMesher.mesh().extractBoundary(1E-12);

    CStlUtils<decimal> aStl(aSurfaceMesh);

    EXPECT_EQ(48, aStl.stlFile().nbFacets());

    aStl.generateConnectivity(1E-12);

    aStl.checkEdges(30.0*3.14 / 180.0);

    aStl.splitFacets(0.5, -1.0);

    aStl.relaxVertices();

    aStl.flipEdges(5.0*3.14 / 180.0);

    EXPECT_EQ(96, aStl.stlFile().nbFacets());

    aStl.save("box_remesh.stl");

}

TEST_F(CTestStlUtils, collapseEdges)
{

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.5, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex6(1.0, 2.0, 0.0);
    CGeoCoordinate<decimal> aVertex7(0.0, 2.0, 0.0);

    CStlFacet<decimal> aFacet1;
    aFacet1.addVertex(aVertex1);
    aFacet1.addVertex(aVertex2);
    aFacet1.addVertex(aVertex4);

    CStlFacet<decimal> aFacet2;
    aFacet2.addVertex(aVertex2);
    aFacet2.addVertex(aVertex3);
    aFacet2.addVertex(aVertex4);

    CStlFacet<decimal> aFacet3;
    aFacet3.addVertex(aVertex1);
    aFacet3.addVertex(aVertex4);
    aFacet3.addVertex(aVertex5);

    CStlFacet<decimal> aFacet4;
    aFacet4.addVertex(aVertex4);
    aFacet4.addVertex(aVertex3);
    aFacet4.addVertex(aVertex6);

    CStlFacet<decimal> aFacet5;
    aFacet5.addVertex(aVertex4);
    aFacet5.addVertex(aVertex6);
    aFacet5.addVertex(aVertex7);

    CStlFacet<decimal> aFacet6;
    aFacet6.addVertex(aVertex5);
    aFacet6.addVertex(aVertex4);
    aFacet6.addVertex(aVertex7);

    CStlUtils<decimal> aStl;

    aStl.addFacet(0, aFacet1);
    aStl.addFacet(1, aFacet2);
    aStl.addFacet(2, aFacet3);
    aStl.addFacet(3, aFacet4);
    aStl.addFacet(4, aFacet5);
    aStl.addFacet(5, aFacet6);

    EXPECT_EQ(6, aStl.stlFile().nbFacets());

    aStl.generateConnectivity(1E-12);
    aStl.collapseEdges(1.0, 5.0*3.14 / 180.0);
    aStl.removeInvalidFacets(1E-12);

    EXPECT_EQ(4, aStl.stlFile().nbFacets());

    aStl.save("plane_collapse.stl");

}

TEST_F(CTestStlUtils, setOrientation1)
{

    CStlUtils<decimal> aStl;

    CGeoCoordinate<decimal> aVertex0(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex1(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(2.0, 1.0, 0.0);

    CStlFacet<decimal> aFacet0;

    aFacet0.addVertex(aVertex0);
    aFacet0.addVertex(aVertex1);
    aFacet0.addVertex(aVertex2);

    aStl.addFacet(0, aFacet0);

    CStlFacet<decimal> aFacet1;

    aFacet1.addVertex(aVertex0);
    aFacet1.addVertex(aVertex2);
    aFacet1.addVertex(aVertex3);

    aStl.addFacet(1, aFacet1);

    CStlFacet<decimal> aFacet2;

    aFacet2.addVertex(aVertex1);
    aFacet2.addVertex(aVertex2);
    aFacet2.addVertex(aVertex5);

    aStl.addFacet(2, aFacet2);

    CStlFacet<decimal> aFacet3;

    aFacet3.addVertex(aVertex1);
    aFacet3.addVertex(aVertex4);
    aFacet3.addVertex(aVertex5);

    aStl.addFacet(3, aFacet3);

    aStl.generateConnectivity(1E-12);

    aStl.setOrientation(0, 1E-12);

    EXPECT_EQ(-1, aStl.stlFile().facet(0).edge(0).neighbor());
    EXPECT_EQ(2, aStl.stlFile().facet(0).edge(1).neighbor());
    EXPECT_EQ(1, aStl.stlFile().facet(0).edge(2).neighbor());

    EXPECT_EQ(-1, aStl.stlFile().facet(0).edge(0).whichVertexNot());
    EXPECT_EQ(2, aStl.stlFile().facet(0).edge(1).whichVertexNot());
    EXPECT_EQ(2, aStl.stlFile().facet(0).edge(2).whichVertexNot());

    EXPECT_EQ(0, aStl.stlFile().facet(2).edge(0).neighbor());
    EXPECT_EQ(3, aStl.stlFile().facet(2).edge(1).neighbor());
    EXPECT_EQ(-1, aStl.stlFile().facet(2).edge(2).neighbor());

    EXPECT_EQ(0, aStl.stlFile().facet(2).edge(0).whichVertexNot());
    EXPECT_EQ(1, aStl.stlFile().facet(2).edge(1).whichVertexNot());
    EXPECT_EQ(-1, aStl.stlFile().facet(2).edge(2).whichVertexNot());

    EXPECT_EQ(-1, aStl.stlFile().facet(3).edge(0).neighbor());
    EXPECT_EQ(-1, aStl.stlFile().facet(3).edge(1).neighbor());
    EXPECT_EQ(2, aStl.stlFile().facet(3).edge(2).neighbor());

    EXPECT_EQ(-1, aStl.stlFile().facet(3).edge(0).whichVertexNot());
    EXPECT_EQ(-1, aStl.stlFile().facet(3).edge(1).whichVertexNot());
    EXPECT_EQ(0, aStl.stlFile().facet(3).edge(2).whichVertexNot());

    aStl.save("plane1_fixnorm.stl");

}

TEST_F(CTestStlUtils, setOrientation2)
{

    CStlUtils<decimal> aStl;

    CGeoCoordinate<decimal> aVertex0(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex1(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(2.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(2.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex6(0.0, 2.0, 0.0);
    CGeoCoordinate<decimal> aVertex7(1.0, 2.0, 0.0);
    CGeoCoordinate<decimal> aVertex8(2.0, 2.0, 0.0);

    CStlFacet<decimal> aFacet0;

    aFacet0.addVertex(aVertex0);
    aFacet0.addVertex(aVertex1);
    aFacet0.addVertex(aVertex2);

    aStl.addFacet(0, aFacet0);

    CStlFacet<decimal> aFacet1;

    aFacet1.addVertex(aVertex2);
    aFacet1.addVertex(aVertex0);
    aFacet1.addVertex(aVertex3);

    aStl.addFacet(1, aFacet1);

    CStlFacet<decimal> aFacet2;

    aFacet2.addVertex(aVertex1);
    aFacet2.addVertex(aVertex2);
    aFacet2.addVertex(aVertex5);

    aStl.addFacet(2, aFacet2);

    CStlFacet<decimal> aFacet3;

    aFacet3.addVertex(aVertex1);
    aFacet3.addVertex(aVertex4);
    aFacet3.addVertex(aVertex5);

    aStl.addFacet(3, aFacet3);

    CStlFacet<decimal> aFacet4;

    aFacet4.addVertex(aVertex3);
    aFacet4.addVertex(aVertex7);
    aFacet4.addVertex(aVertex6);

    aStl.addFacet(4, aFacet4);

    CStlFacet<decimal> aFacet5;

    aFacet5.addVertex(aVertex2);
    aFacet5.addVertex(aVertex3);
    aFacet5.addVertex(aVertex7);

    aStl.addFacet(5, aFacet5);

    CStlFacet<decimal> aFacet6;

    aFacet6.addVertex(aVertex2);
    aFacet6.addVertex(aVertex7);
    aFacet6.addVertex(aVertex8);

    aStl.addFacet(6, aFacet6);

    CStlFacet<decimal> aFacet7;

    aFacet7.addVertex(aVertex2);
    aFacet7.addVertex(aVertex5);
    aFacet7.addVertex(aVertex8);

    aStl.addFacet(7, aFacet7);

    aStl.generateConnectivity(1E-12);

    aStl.setOrientation(0, 1E-12);

    EXPECT_EQ(0, aStl.stlFile().facet(2).edge(0).neighbor());
    EXPECT_EQ(3, aStl.stlFile().facet(2).edge(1).neighbor());
    EXPECT_EQ(7, aStl.stlFile().facet(2).edge(2).neighbor());

    EXPECT_EQ(0, aStl.stlFile().facet(2).edge(0).whichVertexNot());
    EXPECT_EQ(1, aStl.stlFile().facet(2).edge(1).whichVertexNot());
    EXPECT_EQ(2, aStl.stlFile().facet(2).edge(2).whichVertexNot());

    EXPECT_EQ(1, aStl.stlFile().facet(5).edge(0).neighbor());
    EXPECT_EQ(6, aStl.stlFile().facet(5).edge(1).neighbor());
    EXPECT_EQ(4, aStl.stlFile().facet(5).edge(2).neighbor());

    EXPECT_EQ(0, aStl.stlFile().facet(5).edge(0).whichVertexNot());
    EXPECT_EQ(0, aStl.stlFile().facet(5).edge(1).whichVertexNot());
    EXPECT_EQ(2, aStl.stlFile().facet(5).edge(2).whichVertexNot());

    aStl.save("plane2_fixnorm.stl");

}

TEST_F(CTestStlUtils, setOrientation3)
{

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

    const Integer nu = 3;
    const Integer nv = 3;

    aBasicMesher.generate(aQuadrilateral, nu, nv, true);

    aBasicMesher.mesh().element(12).invert();
    aBasicMesher.mesh().element(13).invert();

    CPosGmsh<decimal> aPosGmsh;
    CPdeField<decimal> aField;

    aField.setMesh(aBasicMesher.mesh());
    aPosGmsh.save(aField, "plane3_fixnorm.msh", "");

    CStlUtils<decimal> aStl(aBasicMesher.mesh());

    aStl.generateConnectivity(1E-12);

    aStl.save("plane3.stl");

    aStl.setOrientation(0, 1E-12);

    EXPECT_EQ(15, aStl.stlFile().facet(12).edge(0).neighbor());
    EXPECT_EQ(13, aStl.stlFile().facet(12).edge(1).neighbor());
    EXPECT_EQ(7, aStl.stlFile().facet(12).edge(2).neighbor());

    EXPECT_EQ(1, aStl.stlFile().facet(12).edge(0).whichVertexNot());
    EXPECT_EQ(1, aStl.stlFile().facet(12).edge(1).whichVertexNot());
    EXPECT_EQ(0, aStl.stlFile().facet(12).edge(2).whichVertexNot());

    EXPECT_EQ(-1, aStl.stlFile().facet(13).edge(0).neighbor());
    EXPECT_EQ(-1, aStl.stlFile().facet(13).edge(1).neighbor());
    EXPECT_EQ(12, aStl.stlFile().facet(13).edge(2).neighbor());

    EXPECT_EQ(-1, aStl.stlFile().facet(13).edge(0).whichVertexNot());
    EXPECT_EQ(-1, aStl.stlFile().facet(13).edge(1).whichVertexNot());
    EXPECT_EQ(0, aStl.stlFile().facet(13).edge(2).whichVertexNot());

    aStl.save("plane3_fixnorm.stl");

}
