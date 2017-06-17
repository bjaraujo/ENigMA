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

#include "FvmMesh.hpp"
#include "FvmMeshSearch.hpp"
#include "MshBasicMesher.hpp"

using namespace ENigMA::fvm;

class CTestFvmMesh : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestFvmMesh, calcVolume1) {

    /*
    1 0.0206083 -0.00442002 0.0089
    2 0.02147799300482443 -0.00452017265334567 0.008806528084091029
    3 0.0214445 -0.004 0.0079416
    4 0.0214445 -0.004 0.0089
    */

    CGeoCoordinate<decimal> aVertex1(0.0206083, -0.00442002, 0.0089);
    CGeoCoordinate<decimal> aVertex2(0.02147799300482443, -0.00452017265334567, 0.008806528084091029);
    CGeoCoordinate<decimal> aVertex3(0.0214445, -0.004, 0.0079416);
    CGeoCoordinate<decimal> aVertex4(0.0214445, -0.004, 0.0089);

    CGeoTetrahedron<decimal> aTetrahedron;

    // Inverted 1 and 2
    aTetrahedron.addVertex(aVertex2);
    aTetrahedron.addVertex(aVertex1);
    aTetrahedron.addVertex(aVertex3);
    aTetrahedron.addVertex(aVertex4);

    aTetrahedron.calculateVolume(true);

    decimal volume1 = aTetrahedron.volume();

    EXPECT_GT(volume1, 0.0);

    CMshMesh<decimal> aMesh;

    CMshNode<decimal> aNode1;
    aNode1 << aVertex1;
    aMesh.addNode(1, aNode1);

    CMshNode<decimal> aNode2;
    aNode2 << aVertex2;
    aMesh.addNode(2, aNode2);

    CMshNode<decimal> aNode3;
    aNode3 << aVertex3;
    aMesh.addNode(3, aNode3);

    CMshNode<decimal> aNode4;
    aNode4 << aVertex4;
    aMesh.addNode(4, aNode4);

    CMshElement<decimal> anElement1(ET_TETRAHEDRON);
    anElement1.addNodeId(2);
    anElement1.addNodeId(1);
    anElement1.addNodeId(3);
    anElement1.addNodeId(4);
    aMesh.addElement(1, anElement1);

    aMesh.generateFaces(1E-6);

    CFvmMesh<decimal> aFvmMesh(aMesh);

    aFvmMesh.controlVolume(1).calculateVolume(true);

    decimal volume2 = aFvmMesh.controlVolume(1).volume();

    EXPECT_GT(volume2, 0.0);

}

TEST_F(CTestFvmMesh, clip1) {

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

    aHexahedron.calculateVolume();

    EXPECT_NEAR(1.0, aHexahedron.volume(), 1E-12);

    CMshBasicMesher<decimal> aBasicMesher;

    const Integer nu = 3;
    const Integer nv = 2;
    const Integer nw = 2;

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    aBasicMesher.mesh().generateFaces(1E-12);

    Integer nbTotalFaces = nu * nv * nw * 6;                                   // Total number of faces
    Integer nbExteriorFaces = nu * nv * 2 + nu * nv * 2 + nv * nw * 2;         // Number of exterior faces

    EXPECT_EQ(nbTotalFaces, aBasicMesher.mesh().nbFaces());
    EXPECT_EQ(nbExteriorFaces, aBasicMesher.mesh().nbBoundaryFaces());

    CFvmMesh<decimal> aFvmMesh(aBasicMesher.mesh());

    Integer nbDuplicateFaces = (nbTotalFaces - nbExteriorFaces);               // Number of duplicate faces (interior)
    Integer nbNonDuplicateInteriorFaces = nbDuplicateFaces / 2;                // Number of non-duplicate faces (interior)
    Integer nbUniqueFaces = nbNonDuplicateInteriorFaces + nbExteriorFaces;

    EXPECT_EQ(nbUniqueFaces, aFvmMesh.nbFaces());

    Integer aControlId = 0;

    CFvmControlVolume<decimal> aControlVolume = aFvmMesh.controlVolume(aControlId);
    
    aControlVolume.calculateVolume();

    EXPECT_NEAR(1.0 / (nu * nv * nw), aControlVolume.volume(), 1E-6);

    aControlVolume.calculateSurfaceArea();

    EXPECT_NEAR(1.0 / (nu * nv) * 2 + 1.0 / (nv * nw) * 2 + 1.0 / (nw * nu) * 2, aControlVolume.surfaceArea(), 1E-6);

    EXPECT_EQ(6, aControlVolume.nbFaces());

    Integer aFaceId;

    CFvmMeshSearch<decimal> aFvmMeshSearch(aFvmMesh);

    aFvmMeshSearch.build();

    CGeoCoordinate<decimal> aCoordinate(0.5, 0.75, 1.0);

    aFvmMeshSearch.findClosestBoundaryFace(aCoordinate, aFaceId, 1E-2);

    EXPECT_EQ(46, aFaceId);

    // Current CV
    aControlId = 5;

    aControlVolume = aFvmMesh.controlVolume(aControlId);

    for (Integer i = 0; i < aControlVolume.nbFaces(); ++i)
    {

        Integer aFaceId = aControlVolume.faceId(i);

        CFvmFace<decimal> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateNormal();

        Integer thisControlVolumeId;
        Integer thisNeighborId;

        if (aFace.hasPair())
        {

            thisNeighborId = aFace.neighborId(aControlId);
            thisControlVolumeId = aFace.controlVolumeId(aControlId);

        }
        else
        {
            thisControlVolumeId = aFace.controlVolumeId(aControlId);
            thisNeighborId = -1;
        }

        EXPECT_EQ(aControlId, thisControlVolumeId);

        if (!aFace.hasPair())
            EXPECT_EQ(-1, thisNeighborId);
        else
            EXPECT_NE(-1, thisNeighborId);

    }

    // Clip current CV
    decimal volumeFactionAct;
    Integer nIterations;

    CGeoNormal<decimal> aNormal(1, 0, 0);

    CFvmFace<decimal> aNewFace;

    // Get highest face id
    Integer aNewFaceId = aFvmMesh.nextFaceId();

    aControlVolume.setClippedFaceId(aNewFaceId);
    aControlVolume.clip(aNormal, 0.5, volumeFactionAct, nIterations, 50, 1E-3);

    // Add new face to mesh
    aNewFace.setControlVolumeId(aControlId);
    aFvmMesh.addFace(aNewFaceId, aNewFace);

    for (Integer i = 0; i < aControlVolume.nbFaces(); ++i)
    {

        Integer aFaceId = aControlVolume.faceId(i);

        CFvmFace<decimal> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateNormal();

        Integer thisControlVolumeId;
        Integer thisNeighborId;

        if (aFace.hasPair())
        {

            thisNeighborId = aFace.neighborId(aControlId);
            thisControlVolumeId = aFace.controlVolumeId(aControlId);

        }
        else
        {
            thisControlVolumeId = aFace.controlVolumeId(aControlId);
            thisNeighborId = -1;
        }

        EXPECT_EQ(aControlId, thisControlVolumeId);

        if (!aFace.hasPair())
            EXPECT_EQ(-1, thisNeighborId);
        else
            EXPECT_NE(-1, thisNeighborId);

    }

}

TEST_F(CTestFvmMesh, clip2) {

    decimal length = 0.02;
    decimal width = 0.01;
    decimal height = 0.001;

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(length, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex3(length, width, 0.0);
    CGeoCoordinate<decimal> aVertex4(0.0, width, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 0.0, height);
    CGeoCoordinate<decimal> aVertex6(length, 0.0, height);
    CGeoCoordinate<decimal> aVertex7(length, width, height);
    CGeoCoordinate<decimal> aVertex8(0.0, width, height);

    CGeoHexahedron<decimal> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    const Integer nu = 10;
    const Integer nv = 5;
    const Integer nw = 1;

    CMshBasicMesher<decimal> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    aBasicMesher.mesh().generateFaces(1E-12);

    CFvmMesh<decimal> aFvmMesh(aBasicMesher.mesh());

    // Clip current CV
    Integer aControlId = 25;

    decimal volumeFactionAct;
    Integer nIterations;

    CGeoNormal<decimal> aNormal(1, 0, 0);

    CFvmControlVolume<decimal> aControlVolume;

    aControlVolume = aFvmMesh.controlVolume(aControlId);

    CFvmFace<decimal> aNewFace;

    // Get highest face id
    Integer aNewFaceId = aFvmMesh.nextFaceId();

    aControlVolume.setClippedFaceId(aNewFaceId);
    aControlVolume.clip(aNormal, 1E-7, volumeFactionAct, nIterations, 50, 1E-12);

    EXPECT_NEAR(1E-7, volumeFactionAct, 1E-5);

    for (Integer i = 0; i < aControlVolume.nbFaces(); ++i)
    {

        Integer aFaceId = aControlVolume.faceId(i);

        aControlVolume.calculateFaceArea(aFaceId);

        EXPECT_GT(aControlVolume.faceArea(aFaceId), 0.0);
    }

}
