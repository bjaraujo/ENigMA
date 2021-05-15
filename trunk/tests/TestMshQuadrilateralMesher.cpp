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

#include "MshQuadrilateralMesher.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::pde;

class CTestMshQuadrilateralMesher : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestMshQuadrilateralMesher, mesh1) {

    CMshMesh<decimal> anEdgeMesh;

    CMshNode<decimal> aNode1;
    aNode1 << 0.0, 0.0, 0.0;
    anEdgeMesh.addNode(0, aNode1);

    CMshNode<decimal> aNode2;
    aNode2 << 1.0, 0.0, 0.0;
    anEdgeMesh.addNode(1, aNode2);

    CMshNode<decimal> aNode3;
    aNode3 << 1.0, 1.0, 0.0;
    anEdgeMesh.addNode(2, aNode3);

    CMshNode<decimal> aNode4;
    aNode4 << 0.0, 1.0, 0.0;
    anEdgeMesh.addNode(3, aNode4);

    CMshElement<decimal> anElement1(ET_BEAM);
    anElement1.addNodeId(0);
    anElement1.addNodeId(1);
    anEdgeMesh.addElement(0, anElement1);

    CMshElement<decimal> anElement2(ET_BEAM);
    anElement2.addNodeId(1);
    anElement2.addNodeId(2);
    anEdgeMesh.addElement(1, anElement2);

    CMshElement<decimal> anElement3(ET_BEAM);
    anElement3.addNodeId(2);
    anElement3.addNodeId(3);
    anEdgeMesh.addElement(2, anElement3);

    CMshElement<decimal> anElement4(ET_BEAM);
    anElement4.addNodeId(3);
    anElement4.addNodeId(0);
    anEdgeMesh.addElement(3, anElement4);

    EXPECT_EQ(4, anEdgeMesh.nbElements());

    CMshQuadrilateralMesher<decimal> aQuadrilateralMesher;

    anEdgeMesh.generateFaces(1E-3);
    aQuadrilateralMesher.remesh(anEdgeMesh, 0.1);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "quad_edge1.msh", "beams");

    aQuadrilateralMesher.generate(anEdgeMesh, 999, 0.1, 0.01, 1.0, 1E-3);

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aQuadrilateralMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "quad_surface1.msh", "quads");

    EXPECT_EQ(100, aSurfaceMesh.nbElements());

}

TEST_F(CTestMshQuadrilateralMesher, mesh2) {

    CMshMesh<decimal> anEdgeMesh;

    CMshNode<decimal> aNode1;
    aNode1 << 0.0, 0.0, 0.0;
    anEdgeMesh.addNode(0, aNode1);

    CMshNode<decimal> aNode2;
    aNode2 << 1.0, 0.0, 0.0;
    anEdgeMesh.addNode(1, aNode2);

    CMshNode<decimal> aNode3;
    aNode3 << 1.0, 1.0, 0.0;
    anEdgeMesh.addNode(2, aNode3);

    CMshNode<decimal> aNode4;
    aNode4 << 0.0, 1.0, 0.0;
    anEdgeMesh.addNode(3, aNode4);

    CMshNode<decimal> aNode5;
    aNode5 << 0.3, 0.3, 0.0;
    anEdgeMesh.addNode(4, aNode5);

    CMshNode<decimal> aNode6;
    aNode6 << 0.3, 0.6, 0.0;
    anEdgeMesh.addNode(5, aNode6);

    CMshNode<decimal> aNode7;
    aNode7 << 0.6, 0.6, 0.0;
    anEdgeMesh.addNode(6, aNode7);

    CMshNode<decimal> aNode8;
    aNode8 << 0.6, 0.3, 0.0;
    anEdgeMesh.addNode(7, aNode8);

    CMshElement<decimal> anElement1(ET_BEAM);
    anElement1.addNodeId(0);
    anElement1.addNodeId(1);
    anEdgeMesh.addElement(0, anElement1);

    CMshElement<decimal> anElement2(ET_BEAM);
    anElement2.addNodeId(1);
    anElement2.addNodeId(2);
    anEdgeMesh.addElement(1, anElement2);

    CMshElement<decimal> anElement3(ET_BEAM);
    anElement3.addNodeId(2);
    anElement3.addNodeId(3);
    anEdgeMesh.addElement(2, anElement3);

    CMshElement<decimal> anElement4(ET_BEAM);
    anElement4.addNodeId(3);
    anElement4.addNodeId(0);
    anEdgeMesh.addElement(3, anElement4);

    CMshElement<decimal> anElement5(ET_BEAM);
    anElement5.addNodeId(4);
    anElement5.addNodeId(5);
    anEdgeMesh.addElement(4, anElement5);

    CMshElement<decimal> anElement6(ET_BEAM);
    anElement6.addNodeId(5);
    anElement6.addNodeId(6);
    anEdgeMesh.addElement(5, anElement6);

    CMshElement<decimal> anElement7(ET_BEAM);
    anElement7.addNodeId(6);
    anElement7.addNodeId(7);
    anEdgeMesh.addElement(6, anElement7);

    CMshElement<decimal> anElement8(ET_BEAM);
    anElement8.addNodeId(7);
    anElement8.addNodeId(4);
    anEdgeMesh.addElement(7, anElement8);

    EXPECT_EQ(8, anEdgeMesh.nbElements());

    CMshQuadrilateralMesher<decimal> aQuadrilateralMesher;

    anEdgeMesh.generateFaces(1E-3);
    aQuadrilateralMesher.remesh(anEdgeMesh, 0.1);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "quad_edge2.msh", "beams");

    aQuadrilateralMesher.generate(anEdgeMesh, 999, 0.1, 0.01, 1.0, 1E-3);

    EXPECT_EQ(52, anEdgeMesh.nbElements());

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aQuadrilateralMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "quad_surface2.msh", "quads");

    EXPECT_EQ(91, aSurfaceMesh.nbElements());

}

TEST_F(CTestMshQuadrilateralMesher, mesh3)
{

    CMshMesh<decimal> anEdgeMesh;

    CMshNode<decimal> aNode1;
    aNode1 << 3.004417696365175, 1.91304347826087, 0;
    anEdgeMesh.addNode(0, aNode1);

    CMshNode<decimal> aNode2;
    aNode2 << 2.628865484319528, 1.91304347826087, 0;
    anEdgeMesh.addNode(1, aNode2);

    CMshNode<decimal> aNode3;
    aNode3 << 2.253313272273881, 1.91304347826087, 0;
    anEdgeMesh.addNode(2, aNode3);

    CMshNode<decimal> aNode4;
    aNode4 << 1.877761060228234, 1.91304347826087, 0;
    anEdgeMesh.addNode(3, aNode4);

    CMshNode<decimal> aNode5;
    aNode5 << 1.502208848182587, 1.91304347826087, 0;
    anEdgeMesh.addNode(4, aNode5);

    CMshNode<decimal> aNode6;
    aNode6 << 1.126656636136941, 1.91304347826087, 0;
    anEdgeMesh.addNode(5, aNode6);

    CMshNode<decimal> aNode7;
    aNode7 << 0.7511044240912936, 1.91304347826087, 0;
    anEdgeMesh.addNode(6, aNode7);

    CMshNode<decimal> aNode8;
    aNode8 << 0.3755522120456468, 1.91304347826087, 0;
    anEdgeMesh.addNode(7, aNode8);

    CMshNode<decimal> aNode9;
    aNode9 << 0, 1.91304347826087, 0;
    anEdgeMesh.addNode(8, aNode9);

    CMshNode<decimal> aNode10;
    aNode10 << 3.004417696365175, 0, 0;
    anEdgeMesh.addNode(9, aNode10);

    CMshNode<decimal> aNode11;
    aNode11 << 3.004417696365175, 0.3826086956521739, 0;
    anEdgeMesh.addNode(10, aNode11);

    CMshNode<decimal> aNode12;
    aNode12 << 3.004417696365175, 0.7652173913043478, 0;
    anEdgeMesh.addNode(11, aNode12);

    CMshNode<decimal> aNode13;
    aNode13 << 3.004417696365175, 1.147826086956522, 0;
    anEdgeMesh.addNode(12, aNode13);

    CMshNode<decimal> aNode14;
    aNode14 << 3.004417696365175, 1.530434782608696, 0;
    anEdgeMesh.addNode(13, aNode14);

    CMshNode<decimal> aNode15;
    aNode15 << 0, 0, 0;
    anEdgeMesh.addNode(14, aNode15);

    CMshNode<decimal> aNode16;
    aNode16 << 0.3755522120456473, 0, 0;
    anEdgeMesh.addNode(15, aNode16);

    CMshNode<decimal> aNode17;
    aNode17 << 0.7511044240912945, 0, 0;
    anEdgeMesh.addNode(16, aNode17);

    CMshNode<decimal> aNode18;
    aNode18 << 1.126656636136942, 0, 0;
    anEdgeMesh.addNode(17, aNode18);

    CMshNode<decimal> aNode19;
    aNode19 << 1.502208848182589, 0, 0;
    anEdgeMesh.addNode(18, aNode19);

    CMshNode<decimal> aNode20;
    aNode20 << 1.877761060228236, 0, 0;
    anEdgeMesh.addNode(19, aNode20);

    CMshNode<decimal> aNode21;
    aNode21 << 2.253313272273883, 0, 0;
    anEdgeMesh.addNode(20, aNode21);

    CMshNode<decimal> aNode22;
    aNode22 << 2.628865484319531, 0, 0;
    anEdgeMesh.addNode(21, aNode22);

    CMshNode<decimal> aNode23;
    aNode23 << 0, 1.530434782608696, 0;
    anEdgeMesh.addNode(22, aNode23);

    CMshNode<decimal> aNode24;
    aNode24 << 0, 1.147826086956522, 0;
    anEdgeMesh.addNode(23, aNode24);

    CMshNode<decimal> aNode25;
    aNode25 << 0, 0.7652173913043476, 0;
    anEdgeMesh.addNode(24, aNode25);

    CMshNode<decimal> aNode26;
    aNode26 << 0, 0.3826086956521738, 0;
    anEdgeMesh.addNode(25, aNode26);

    EXPECT_EQ(26, anEdgeMesh.nbNodes());

    CMshElement<decimal> anElement1(ET_BEAM);
    anElement1.addNodeId(0);
    anElement1.addNodeId(1);
    anEdgeMesh.addElement(0, anElement1);

    CMshElement<decimal> anElement2(ET_BEAM);
    anElement2.addNodeId(1);
    anElement2.addNodeId(2);
    anEdgeMesh.addElement(1, anElement2);

    CMshElement<decimal> anElement3(ET_BEAM);
    anElement3.addNodeId(2);
    anElement3.addNodeId(3);
    anEdgeMesh.addElement(2, anElement3);

    CMshElement<decimal> anElement4(ET_BEAM);
    anElement4.addNodeId(3);
    anElement4.addNodeId(4);
    anEdgeMesh.addElement(3, anElement4);

    CMshElement<decimal> anElement5(ET_BEAM);
    anElement5.addNodeId(4);
    anElement5.addNodeId(5);
    anEdgeMesh.addElement(4, anElement5);

    CMshElement<decimal> anElement6(ET_BEAM);
    anElement6.addNodeId(5);
    anElement6.addNodeId(6);
    anEdgeMesh.addElement(5, anElement6);

    CMshElement<decimal> anElement7(ET_BEAM);
    anElement7.addNodeId(6);
    anElement7.addNodeId(7);
    anEdgeMesh.addElement(6, anElement7);

    CMshElement<decimal> anElement8(ET_BEAM);
    anElement8.addNodeId(7);
    anElement8.addNodeId(8);
    anEdgeMesh.addElement(7, anElement8);

    CMshElement<decimal> anElement9(ET_BEAM);
    anElement9.addNodeId(9);
    anElement9.addNodeId(10);
    anEdgeMesh.addElement(8, anElement9);

    CMshElement<decimal> anElement10(ET_BEAM);
    anElement10.addNodeId(10);
    anElement10.addNodeId(11);
    anEdgeMesh.addElement(9, anElement10);

    CMshElement<decimal> anElement11(ET_BEAM);
    anElement11.addNodeId(11);
    anElement11.addNodeId(12);
    anEdgeMesh.addElement(10, anElement11);

    CMshElement<decimal> anElement12(ET_BEAM);
    anElement12.addNodeId(12);
    anElement12.addNodeId(13);
    anEdgeMesh.addElement(11, anElement12);

    CMshElement<decimal> anElement13(ET_BEAM);
    anElement13.addNodeId(13);
    anElement13.addNodeId(0);
    anEdgeMesh.addElement(12, anElement13);

    CMshElement<decimal> anElement14(ET_BEAM);
    anElement14.addNodeId(14);
    anElement14.addNodeId(15);
    anEdgeMesh.addElement(13, anElement14);

    CMshElement<decimal> anElement15(ET_BEAM);
    anElement15.addNodeId(15);
    anElement15.addNodeId(16);
    anEdgeMesh.addElement(14, anElement15);

    CMshElement<decimal> anElement16(ET_BEAM);
    anElement16.addNodeId(16);
    anElement16.addNodeId(17);
    anEdgeMesh.addElement(15, anElement16);

    CMshElement<decimal> anElement17(ET_BEAM);
    anElement17.addNodeId(17);
    anElement17.addNodeId(18);
    anEdgeMesh.addElement(16, anElement17);

    CMshElement<decimal> anElement18(ET_BEAM);
    anElement18.addNodeId(18);
    anElement18.addNodeId(19);
    anEdgeMesh.addElement(17, anElement18);

    CMshElement<decimal> anElement19(ET_BEAM);
    anElement19.addNodeId(19);
    anElement19.addNodeId(20);
    anEdgeMesh.addElement(18, anElement19);

    CMshElement<decimal> anElement20(ET_BEAM);
    anElement20.addNodeId(20);
    anElement20.addNodeId(21);
    anEdgeMesh.addElement(19, anElement20);

    CMshElement<decimal> anElement21(ET_BEAM);
    anElement21.addNodeId(21);
    anElement21.addNodeId(9);
    anEdgeMesh.addElement(20, anElement21);

    CMshElement<decimal> anElement22(ET_BEAM);
    anElement22.addNodeId(8);
    anElement22.addNodeId(22);
    anEdgeMesh.addElement(21, anElement22);

    CMshElement<decimal> anElement23(ET_BEAM);
    anElement23.addNodeId(22);
    anElement23.addNodeId(23);
    anEdgeMesh.addElement(22, anElement23);

    CMshElement<decimal> anElement24(ET_BEAM);
    anElement24.addNodeId(23);
    anElement24.addNodeId(24);
    anEdgeMesh.addElement(23, anElement24);

    CMshElement<decimal> anElement25(ET_BEAM);
    anElement25.addNodeId(24);
    anElement25.addNodeId(25);
    anEdgeMesh.addElement(24, anElement25);

    CMshElement<decimal> anElement26(ET_BEAM);
    anElement26.addNodeId(25);
    anElement26.addNodeId(14);
    anEdgeMesh.addElement(25, anElement26);

    EXPECT_EQ(26, anEdgeMesh.nbElements());

    CMshQuadrilateralMesher<decimal> aQuadrilateralMesher;

    anEdgeMesh.generateFaces(1E-3);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "quad_edge3.msh", "beams");

    aQuadrilateralMesher.generate(anEdgeMesh, 999, 0.5, 0.5, 0.5, 1E-3);

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aQuadrilateralMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "quad_surface3.msh", "quads");

    EXPECT_EQ(35, aSurfaceMesh.nbElements());
}

TEST_F(CTestMshQuadrilateralMesher, flip1)
{

    CMshMesh<decimal> aSurfaceMesh;

    CMshNode<decimal> aNode1;
    aNode1 << 0, 0, 0;
    aSurfaceMesh.addNode(0, aNode1);

    CMshNode<decimal> aNode2;
    aNode2 << 1, 0, 0;
    aSurfaceMesh.addNode(1, aNode2);

    CMshNode<decimal> aNode3;
    aNode3 << 1, 1, 0;
    aSurfaceMesh.addNode(2, aNode3);

    CMshNode<decimal> aNode4;
    aNode4 << 0, 1, 0;
    aSurfaceMesh.addNode(3, aNode4);

    CMshNode<decimal> aNode5;
    aNode5 << 1.05, 0.5, 0;
    aSurfaceMesh.addNode(4, aNode5);

    CMshElement<decimal> anElement1(ET_QUADRILATERAL);
    anElement1.addNodeId(0);
    anElement1.addNodeId(1);
    anElement1.addNodeId(2);
    anElement1.addNodeId(3);
    aSurfaceMesh.addElement(0, anElement1);

    CMshElement<decimal> anElement2(ET_TRIANGLE);
    anElement2.addNodeId(1);
    anElement2.addNodeId(4);
    anElement2.addNodeId(2);
    aSurfaceMesh.addElement(1, anElement2);

    EXPECT_EQ(0, aSurfaceMesh.element(0).nodeId(0));
    EXPECT_EQ(1, aSurfaceMesh.element(0).nodeId(1));
    EXPECT_EQ(2, aSurfaceMesh.element(0).nodeId(2));
    EXPECT_EQ(3, aSurfaceMesh.element(0).nodeId(3));

    CMshQuadrilateralMesher<decimal> aQuadrilateralMesher;

    aSurfaceMesh.generateFaces(1E-3);

    aQuadrilateralMesher.flipEdges(aSurfaceMesh, 1E-3);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "quad_surface4.msh", "quads");

    EXPECT_EQ(2, aSurfaceMesh.element(0).nodeId(0));
    EXPECT_EQ(3, aSurfaceMesh.element(0).nodeId(1));
    EXPECT_EQ(0, aSurfaceMesh.element(0).nodeId(2));
    EXPECT_EQ(4, aSurfaceMesh.element(0).nodeId(3));
}

TEST_F(CTestMshQuadrilateralMesher, flip2)
{

    CMshMesh<decimal> aSurfaceMesh;

    CMshNode<decimal> aNode1;
    aNode1 << 0, 0, 0;
    aSurfaceMesh.addNode(0, aNode1);

    CMshNode<decimal> aNode2;
    aNode2 << 1, 0, 0;
    aSurfaceMesh.addNode(1, aNode2);

    CMshNode<decimal> aNode3;
    aNode3 << 1, 1, 0;
    aSurfaceMesh.addNode(2, aNode3);

    CMshNode<decimal> aNode4;
    aNode4 << 0, 1, 0;
    aSurfaceMesh.addNode(3, aNode4);

    CMshNode<decimal> aNode5;
    aNode5 << 0.5, 1.05, 0;
    aSurfaceMesh.addNode(4, aNode5);

    CMshElement<decimal> anElement1(ET_QUADRILATERAL);
    anElement1.addNodeId(0);
    anElement1.addNodeId(1);
    anElement1.addNodeId(2);
    anElement1.addNodeId(3);
    aSurfaceMesh.addElement(0, anElement1);

    CMshElement<decimal> anElement2(ET_TRIANGLE);
    anElement2.addNodeId(2);
    anElement2.addNodeId(4);
    anElement2.addNodeId(3);
    aSurfaceMesh.addElement(1, anElement2);

    EXPECT_EQ(0, aSurfaceMesh.element(0).nodeId(0));
    EXPECT_EQ(1, aSurfaceMesh.element(0).nodeId(1));
    EXPECT_EQ(2, aSurfaceMesh.element(0).nodeId(2));
    EXPECT_EQ(3, aSurfaceMesh.element(0).nodeId(3));

    CMshQuadrilateralMesher<decimal> aQuadrilateralMesher;

    aSurfaceMesh.generateFaces(1E-3);
    aQuadrilateralMesher.flipEdges(aSurfaceMesh, 1E-3);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(aQuadrilateralMesher.mesh());
    aPosGmsh.save(T, "quad_surface5.msh", "quads");

    EXPECT_EQ(3, aSurfaceMesh.element(0).nodeId(0));
    EXPECT_EQ(0, aSurfaceMesh.element(0).nodeId(1));
    EXPECT_EQ(1, aSurfaceMesh.element(0).nodeId(2));
    EXPECT_EQ(4, aSurfaceMesh.element(0).nodeId(3));
}
