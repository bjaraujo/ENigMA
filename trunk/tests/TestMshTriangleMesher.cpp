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

#include "MshTriangleMesher.hpp"
#include "PdeField.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::pde;

class CTestMshTriangleMesher : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestMshTriangleMesher, mesh1) {

    CMshMesh<decimal> anEdgeMesh;

    CMshNode<decimal> aNode1;
    CMshNode<decimal> aNode2;
    CMshNode<decimal> aNode3;
    CMshNode<decimal> aNode4;
    CMshNode<decimal> aNode5;
    CMshNode<decimal> aNode6;
    CMshNode<decimal> aNode7;
    CMshNode<decimal> aNode8;

    aNode1 << 2.0270543283699989, 1.3490401401130008, 0;
    aNode2 << 2.0270543283699989, 2.6634868158110017, 0;
    aNode3 << 1.3698309905200003, 2.0062634779620012, 0;
    aNode4 << 2.6842776662166674, 2.0062634779620012, 0;
    aNode5 << 2.0270543284199927, 3.8490399419350005, 0;
    aNode6 << 2.0270543251100004, 0.16347718345100049, 0;
    aNode7 << 3.8698309905200001, 2.0062634779619994, 0;
    aNode8 << 0.18427766621999805, 2.0062634779619994, 0;

    anEdgeMesh.addNode(1, aNode1);
    anEdgeMesh.addNode(2, aNode2);
    anEdgeMesh.addNode(3, aNode3);
    anEdgeMesh.addNode(4, aNode4);
    anEdgeMesh.addNode(5, aNode5);
    anEdgeMesh.addNode(6, aNode6);
    anEdgeMesh.addNode(7, aNode7);
    anEdgeMesh.addNode(8, aNode8);

    CMshElement<decimal> anElement1(ET_BEAM);
    CMshElement<decimal> anElement2(ET_BEAM);
    CMshElement<decimal> anElement3(ET_BEAM);
    CMshElement<decimal> anElement4(ET_BEAM);
    CMshElement<decimal> anElement5(ET_BEAM);
    CMshElement<decimal> anElement6(ET_BEAM);
    CMshElement<decimal> anElement7(ET_BEAM);
    CMshElement<decimal> anElement8(ET_BEAM);

    // 1 1 0 1 3
    anElement1.addNodeId(1);
    anElement1.addNodeId(3);

    // 2 1 0 3 2
    anElement2.addNodeId(3);
    anElement2.addNodeId(2);

    // 3 1 0 2 4
    anElement3.addNodeId(2);
    anElement3.addNodeId(4);

    // 4 1 0 4 1
    anElement4.addNodeId(4);
    anElement4.addNodeId(1);

    // 5 1 0 7 5
    anElement5.addNodeId(7);
    anElement5.addNodeId(5);

    // 6 1 0 6 7
    anElement6.addNodeId(6);
    anElement6.addNodeId(7);

    // 7 1 0 8 6
    anElement7.addNodeId(8);
    anElement7.addNodeId(6);

    // 8 1 0 5 8
    anElement8.addNodeId(5);
    anElement8.addNodeId(8);

    anEdgeMesh.addElement(1, anElement1);
    anEdgeMesh.addElement(2, anElement2);
    anEdgeMesh.addElement(3, anElement3);
    anEdgeMesh.addElement(4, anElement4);
    anEdgeMesh.addElement(5, anElement5);
    anEdgeMesh.addElement(6, anElement6);
    anEdgeMesh.addElement(7, anElement7);
    anEdgeMesh.addElement(8, anElement8);

    EXPECT_EQ(8, anEdgeMesh.nbElements());

    CMshTriangleMesher<decimal> aTriangleMesher;

    anEdgeMesh.generateFaces(1E-3);
    aTriangleMesher.remesh(anEdgeMesh, 1.0);
    
    aTriangleMesher.generate(anEdgeMesh, 99, 1.0, 0.1, 1E-3);

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aTriangleMesher.mesh();

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "tri_edge1.msh", "beams");

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tri_surface1.msh", "tris");

    EXPECT_EQ(16, anEdgeMesh.nbElements());
    EXPECT_EQ(16, aSurfaceMesh.nbElements());

}

TEST_F(CTestMshTriangleMesher, mesh2) {

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

    CMshTriangleMesher<decimal> aTriangleMesher;

    anEdgeMesh.generateFaces(1E-3);
    aTriangleMesher.remesh(anEdgeMesh, 0.1);

    aTriangleMesher.flipEdges();
    aTriangleMesher.relaxNodes();

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "tri_edge2.msh", "beams");

    aTriangleMesher.generate(anEdgeMesh, 999, 0.1, 0.1, 1E-3);

    EXPECT_EQ(40, anEdgeMesh.nbElements());

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aTriangleMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tri_surface2.msh", "tris");

    EXPECT_EQ(222, aSurfaceMesh.nbElements());

}

TEST_F(CTestMshTriangleMesher, mesh3) {

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

    CMshTriangleMesher<decimal> aTriangleMesher;

    anEdgeMesh.generateFaces(1E-3);
    aTriangleMesher.remesh(anEdgeMesh, 0.1);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "tri_edge3.msh", "beams");

    aTriangleMesher.generate(anEdgeMesh, 999, 0.1, 0.1, 1E-3);

    aTriangleMesher.flipEdges();
    aTriangleMesher.relaxNodes();

    EXPECT_EQ(52, anEdgeMesh.nbElements());

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aTriangleMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tri_surface3.msh", "tris");

    EXPECT_EQ(216, aSurfaceMesh.nbElements());

}

TEST_F(CTestMshTriangleMesher, mesh4) {

    CMshMesh<decimal> anEdgeMesh;

    CMshNode<decimal> aNode1;
    CMshNode<decimal> aNode2;
    CMshNode<decimal> aNode3;
    CMshNode<decimal> aNode4;
    CMshNode<decimal> aNode5;
    CMshNode<decimal> aNode6;

    aNode1 << 0.0387116, 0, 0;
    aNode2 << 0.0387116, 0.203556, 0;
    aNode3 << 0.0387116, 0.101778, 0;
    aNode4 << 0.000482906, 0.203556, 0;
    aNode5 << -2.09856e-021, -7.06226e-019, 0;
    aNode6 << 0.000250622, 0.101778, 0;

    anEdgeMesh.addNode(1, aNode1);
    anEdgeMesh.addNode(2, aNode2);
    anEdgeMesh.addNode(3, aNode3);
    anEdgeMesh.addNode(4, aNode4);
    anEdgeMesh.addNode(5, aNode5);
    anEdgeMesh.addNode(6, aNode6);

    CMshElement<decimal> anElement1(ET_BEAM);
    CMshElement<decimal> anElement2(ET_BEAM);
    CMshElement<decimal> anElement3(ET_BEAM);
    CMshElement<decimal> anElement4(ET_BEAM);
    CMshElement<decimal> anElement5(ET_BEAM);
    CMshElement<decimal> anElement6(ET_BEAM);

    anElement1.addNodeId(1);
    anElement1.addNodeId(3);
    anEdgeMesh.addElement(1, anElement1);

    anElement2.addNodeId(3);
    anElement2.addNodeId(2);
    anEdgeMesh.addElement(2, anElement2);

    anElement3.addNodeId(2);
    anElement3.addNodeId(4);
    anEdgeMesh.addElement(3, anElement3);

    anElement4.addNodeId(4);
    anElement4.addNodeId(6);
    anEdgeMesh.addElement(4, anElement4);

    anElement5.addNodeId(6);
    anElement5.addNodeId(5);
    anEdgeMesh.addElement(5, anElement5);

    anElement6.addNodeId(5);
    anElement6.addNodeId(1);
    anEdgeMesh.addElement(6, anElement6);

    EXPECT_EQ(6, anEdgeMesh.nbElements());

    CMshTriangleMesher<decimal> aTriangleMesher;

    anEdgeMesh.generateFaces(1E-3);
    aTriangleMesher.remesh(anEdgeMesh, 2.0);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "tri_edge4.msh", "beams");

    aTriangleMesher.generate(anEdgeMesh, 99, 2.0, 0.1, 1E-3);

    aTriangleMesher.flipEdges();
    aTriangleMesher.relaxNodes();

    EXPECT_EQ(6, anEdgeMesh.nbElements());

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aTriangleMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tri_surface4.msh", "tris");

    EXPECT_EQ(4, aSurfaceMesh.nbElements());

}

TEST_F(CTestMshTriangleMesher, mesh5) {

    CMshMesh<decimal> anEdgeMesh;

    CMshNode<decimal> aNode1(0, 0, 0);
    CMshNode<decimal> aNode2(100, 0, 0);
    CMshNode<decimal> aNode3(100, 100, 0);
    CMshNode<decimal> aNode4(0, 100, 0);
    CMshNode<decimal> aNode5(10, 0, 0);
    CMshNode<decimal> aNode6(28, 40, 0);
    CMshNode<decimal> aNode7(28, 0, 0);
    CMshNode<decimal> aNode8(40, 0, 0);
    CMshNode<decimal> aNode9(58, 52, 0);
    CMshNode<decimal> aNode10(58, 0, 0);
    CMshNode<decimal> aNode11(50, 88, 0);
    CMshNode<decimal> aNode12(79, 65, 0);

    anEdgeMesh.addNode(1, aNode1);
    anEdgeMesh.addNode(2, aNode2);
    anEdgeMesh.addNode(3, aNode3);
    anEdgeMesh.addNode(4, aNode4);
    anEdgeMesh.addNode(5, aNode5);
    anEdgeMesh.addNode(6, aNode6);
    anEdgeMesh.addNode(7, aNode7);
    anEdgeMesh.addNode(8, aNode8);
    anEdgeMesh.addNode(9, aNode9);
    anEdgeMesh.addNode(10, aNode10);
    anEdgeMesh.addNode(11, aNode11);
    anEdgeMesh.addNode(12, aNode12);

    CMshElement<decimal> anElement1(ET_BEAM);
    CMshElement<decimal> anElement2(ET_BEAM);
    CMshElement<decimal> anElement3(ET_BEAM);
    CMshElement<decimal> anElement4(ET_BEAM);
    CMshElement<decimal> anElement5(ET_BEAM);
    CMshElement<decimal> anElement6(ET_BEAM);
    CMshElement<decimal> anElement7(ET_BEAM);
    CMshElement<decimal> anElement8(ET_BEAM);
    CMshElement<decimal> anElement9(ET_BEAM);
    CMshElement<decimal> anElement10(ET_BEAM);
    CMshElement<decimal> anElement11(ET_BEAM);
    CMshElement<decimal> anElement12(ET_BEAM);

    anElement1.addNodeId(4);
    anElement1.addNodeId(1);
    anEdgeMesh.addElement(1, anElement1);

    anElement2.addNodeId(1);
    anElement2.addNodeId(5);
    anEdgeMesh.addElement(2, anElement2);

    anElement3.addNodeId(5);
    anElement3.addNodeId(6);
    anEdgeMesh.addElement(3, anElement3);

    anElement4.addNodeId(6);
    anElement4.addNodeId(7);
    anEdgeMesh.addElement(4, anElement4);

    anElement5.addNodeId(7);
    anElement5.addNodeId(8);
    anEdgeMesh.addElement(5, anElement5);

    anElement6.addNodeId(8);
    anElement6.addNodeId(9);
    anEdgeMesh.addElement(6, anElement6);

    anElement7.addNodeId(9);
    anElement7.addNodeId(10);
    anEdgeMesh.addElement(7, anElement7);

    anElement8.addNodeId(10);
    anElement8.addNodeId(2);
    anEdgeMesh.addElement(8, anElement8);

    anElement9.addNodeId(2);
    anElement9.addNodeId(3);
    anEdgeMesh.addElement(9, anElement9);

    anElement10.addNodeId(3);
    anElement10.addNodeId(4);
    anEdgeMesh.addElement(10, anElement10);

    anElement11.addNodeId(11);
    anElement11.addNodeId(12);
    anEdgeMesh.addElement(12, anElement11);

    anElement12.addNodeId(12);
    anElement12.addNodeId(11);
    anEdgeMesh.addElement(13, anElement12);

    EXPECT_EQ(12, anEdgeMesh.nbElements());

    CMshTriangleMesher<decimal> aTriangleMesher;

    anEdgeMesh.generateFaces(1E-3);

    anEdgeMesh.node(5).x() = 20;
    anEdgeMesh.node(8).x() = 50;

    aTriangleMesher.remesh(anEdgeMesh, 10);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "tri_edge5.msh", "beams");

    aTriangleMesher.generate(anEdgeMesh, 999, 10, 0.1, 1E-3);

    aTriangleMesher.flipEdges();
    aTriangleMesher.relaxNodes();

    EXPECT_EQ(64, anEdgeMesh.nbElements());

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aTriangleMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tri_surface5.msh", "tris");

    EXPECT_EQ(aSurfaceMesh.nbElements(), 220);

}

TEST_F(CTestMshTriangleMesher, mesh6) {

    CMshMesh<decimal> anEdgeMesh;

    CMshNode<decimal> aNode1(0, 0, 0);
    CMshNode<decimal> aNode2(100, 0, 0);
    CMshNode<decimal> aNode3(100, 100, 0);
    CMshNode<decimal> aNode4(0, 100, 0);
    CMshNode<decimal> aNode5(28, 0, 0);
    CMshNode<decimal> aNode6(28, 40, 0);
    CMshNode<decimal> aNode7(28, 0, 0);
    CMshNode<decimal> aNode8(58, 0, 0);
    CMshNode<decimal> aNode9(58, 52, 0);
    CMshNode<decimal> aNode10(58, 0, 0);
    CMshNode<decimal> aNode11(50, 88, 0);
    CMshNode<decimal> aNode12(79, 65, 0);

    anEdgeMesh.addNode(1, aNode1);
    anEdgeMesh.addNode(2, aNode2);
    anEdgeMesh.addNode(3, aNode3);
    anEdgeMesh.addNode(4, aNode4);
    anEdgeMesh.addNode(5, aNode5);
    anEdgeMesh.addNode(6, aNode6);
    anEdgeMesh.addNode(7, aNode7);
    anEdgeMesh.addNode(8, aNode8);
    anEdgeMesh.addNode(9, aNode9);
    anEdgeMesh.addNode(10, aNode10);
    anEdgeMesh.addNode(11, aNode11);
    anEdgeMesh.addNode(12, aNode12);

    CMshElement<decimal> anElement1(ET_BEAM);
    CMshElement<decimal> anElement2(ET_BEAM);
    CMshElement<decimal> anElement3(ET_BEAM);
    CMshElement<decimal> anElement4(ET_BEAM);
    CMshElement<decimal> anElement5(ET_BEAM);
    CMshElement<decimal> anElement6(ET_BEAM);
    CMshElement<decimal> anElement7(ET_BEAM);
    CMshElement<decimal> anElement8(ET_BEAM);
    CMshElement<decimal> anElement9(ET_BEAM);
    CMshElement<decimal> anElement10(ET_BEAM);
    CMshElement<decimal> anElement11(ET_BEAM);
    CMshElement<decimal> anElement12(ET_BEAM);

    anElement1.addNodeId(4);
    anElement1.addNodeId(1);
    anEdgeMesh.addElement(1, anElement1);

    anElement2.addNodeId(1);
    anElement2.addNodeId(5);
    anEdgeMesh.addElement(2, anElement2);

    anElement3.addNodeId(5);
    anElement3.addNodeId(8);
    anEdgeMesh.addElement(3, anElement3);

    anElement4.addNodeId(8);
    anElement4.addNodeId(2);
    anEdgeMesh.addElement(4, anElement4);

    anElement5.addNodeId(2);
    anElement5.addNodeId(3);
    anEdgeMesh.addElement(5, anElement5);

    anElement6.addNodeId(3);
    anElement6.addNodeId(4);
    anEdgeMesh.addElement(6, anElement6);

    anElement7.addNodeId(6);
    anElement7.addNodeId(7);
    anEdgeMesh.addElement(7, anElement7);

    anElement8.addNodeId(7);
    anElement8.addNodeId(6);
    anEdgeMesh.addElement(8, anElement8);

    anElement9.addNodeId(9);
    anElement9.addNodeId(10);
    anEdgeMesh.addElement(9, anElement9);

    anElement10.addNodeId(10);
    anElement10.addNodeId(9);
    anEdgeMesh.addElement(10, anElement10);

    anElement11.addNodeId(11);
    anElement11.addNodeId(12);
    anEdgeMesh.addElement(12, anElement11);

    anElement12.addNodeId(12);
    anElement12.addNodeId(11);
    anEdgeMesh.addElement(13, anElement12);

    EXPECT_EQ(12, anEdgeMesh.nbElements());

    CMshTriangleMesher<decimal> aTriangleMesher;

    anEdgeMesh.node(6).z() = 10;
    anEdgeMesh.node(7).z() = 10;

    anEdgeMesh.node(9).z() = 10;
    anEdgeMesh.node(10).z() = 10;

    anEdgeMesh.node(11).z() = 10;
    anEdgeMesh.node(12).z() = 10;

    anEdgeMesh.generateFaces(1E-3);

    anEdgeMesh.node(6).z() = 0;
    anEdgeMesh.node(7).z() = 0;

    anEdgeMesh.node(9).z() = 0;
    anEdgeMesh.node(10).z() = 0;

    anEdgeMesh.node(11).z() = 0;
    anEdgeMesh.node(12).z() = 0;

    aTriangleMesher.remesh(anEdgeMesh, 10);

    anEdgeMesh.mergeNodes(1E-3);

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "tri_edge6.msh", "beams");

    aTriangleMesher.generate(anEdgeMesh, 999, 10, 0.1, 1E-3);
    aTriangleMesher.applyFixedBoundary(anEdgeMesh, 1E-3);

    aTriangleMesher.flipEdges();
    aTriangleMesher.applyFixedBoundary(anEdgeMesh, 1E-3);

    aTriangleMesher.relaxNodes();

    aTriangleMesher.flipEdges();
    aTriangleMesher.applyFixedBoundary(anEdgeMesh, 1E-3);

    EXPECT_EQ(66, anEdgeMesh.nbElements());

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aTriangleMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tri_surface6.msh", "tris");

    EXPECT_EQ(aSurfaceMesh.nbElements(), 238);

}

TEST_F(CTestMshTriangleMesher, mesh7) {

    CMshMesh<decimal> anEdgeMesh;

    CMshNode<decimal> aNode1(0, 0, 0);
    CMshNode<decimal> aNode2(100, 0, 0);
    CMshNode<decimal> aNode3(100, 100, 0);
    CMshNode<decimal> aNode4(0, 100, 0);

    anEdgeMesh.addNode(1, aNode1);
    anEdgeMesh.addNode(2, aNode2);
    anEdgeMesh.addNode(3, aNode3);
    anEdgeMesh.addNode(4, aNode4);

    CMshElement<decimal> anElement1(ET_BEAM);
    CMshElement<decimal> anElement2(ET_BEAM);
    CMshElement<decimal> anElement3(ET_BEAM);
    CMshElement<decimal> anElement4(ET_BEAM);

    anElement1.addNodeId(1);
    anElement1.addNodeId(2);
    anEdgeMesh.addElement(1, anElement1);

    anElement2.addNodeId(2);
    anElement2.addNodeId(3);
    anEdgeMesh.addElement(2, anElement2);

    anElement3.addNodeId(3);
    anElement3.addNodeId(4);
    anEdgeMesh.addElement(3, anElement3);

    anElement4.addNodeId(4);
    anElement4.addNodeId(1);
    anEdgeMesh.addElement(4, anElement4);

    // Inner circle
    decimal r = 5;
    const decimal pi = std::acos(-1.0);

    for (Integer i = 0; i < 16; i++)
    {
        CMshNode<decimal> aNode(60 + r * cos(2 * pi - i * 2 * pi / 16), 60 + r * sin(2 * pi - i * 2 * pi / 16), 0);
        anEdgeMesh.addNode(5 + i, aNode);

        CMshElement<decimal> anElement(ET_BEAM);

        if (i < 15)
        {
            anElement.addNodeId(5 + i);
            anElement.addNodeId(5 + i + 1);
        }
        else
        {
            anElement.addNodeId(5 + i);
            anElement.addNodeId(5 + 0);
        }

        anEdgeMesh.addElement(5 + i, anElement);

    }

    anEdgeMesh.mergeNodes(1E-2);
    anEdgeMesh.generateFaces(1E-2);

    CMshTriangleMesher<decimal> aTriangleMesher;

    aTriangleMesher.remesh(anEdgeMesh, 2.0);
    aTriangleMesher.generate(anEdgeMesh, 9999, 2.0, 0.1, 1E-4);
    
    aTriangleMesher.flipEdges();
    aTriangleMesher.relaxNodes();

    CPdeField<decimal> T;
    CPosGmsh<decimal> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "tri_edge7.msh", "beams");

    CMshMesh<decimal> aSurfaceMesh;
    aSurfaceMesh = aTriangleMesher.mesh();

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "tri_surface7.msh", "tris");

    EXPECT_EQ(5790, aSurfaceMesh.nbElements());

}
