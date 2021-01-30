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

