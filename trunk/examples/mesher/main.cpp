// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iomanip>
#include <iostream>

#include <time.h>

#include "MshTriangleMesher.hpp"
#include "MshQuadrilateralMesher.hpp"
#include "MshTetrahedronMesher.hpp"
#include "PosGmsh.hpp"
#include "StlUtils.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::pde;
using namespace ENigMA::stl;

CMshMesh<float> GenerateBoundary(const float tol)
{
    CMshMesh<float> anEdgeMesh;

    CMshNode<float> aNode1;
    aNode1 << 0.0, 0.0, 0.0;
    anEdgeMesh.addNode(0, aNode1);

    CMshNode<float> aNode2;
    aNode2 << 1.0, 0.0, 0.0;
    anEdgeMesh.addNode(1, aNode2);

    CMshNode<float> aNode3;
    aNode3 << 1.0, 1.0, 0.0;
    anEdgeMesh.addNode(2, aNode3);

    CMshNode<float> aNode4;
    aNode4 << 0.0, 1.0, 0.0;
    anEdgeMesh.addNode(3, aNode4);

    CMshElement<float> anElement1(ET_BEAM);
    anElement1.addNodeId(0);
    anElement1.addNodeId(1);
    anEdgeMesh.addElement(0, anElement1);

    CMshElement<float> anElement2(ET_BEAM);
    anElement2.addNodeId(1);
    anElement2.addNodeId(2);
    anEdgeMesh.addElement(1, anElement2);

    CMshElement<float> anElement3(ET_BEAM);
    anElement3.addNodeId(2);
    anElement3.addNodeId(3);
    anEdgeMesh.addElement(2, anElement3);

    CMshElement<float> anElement4(ET_BEAM);
    anElement4.addNodeId(3);
    anElement4.addNodeId(0);
    anEdgeMesh.addElement(3, anElement4);

    anEdgeMesh.generateFaces(tol);
    anEdgeMesh.mergeNodes(tol);

    return anEdgeMesh;
}

CMshMesh<float> LoadBoundary(std::string fileName, const float tol)
{
    CMshMesh<float> anEdgeMesh;

    CPdeField<float> T;
    CPosGmsh<float> aPosGmsh;

    aPosGmsh.load(T, fileName);
    anEdgeMesh = T.mesh();

    anEdgeMesh.mergeNodes(tol);
    anEdgeMesh.mergeElements(tol);
    anEdgeMesh.renumber();
    anEdgeMesh.generateFaces(tol);

    return anEdgeMesh;
}

void GenerateMesh(CMshMesh<float> anEdgeMesh, const float meshSize, const int maxEle, const float tol)
{
    CMshMesh<float> aSurfaceMesh;

    CPdeField<float> T;
    CPosGmsh<float> aPosGmsh;

    CMshTriangleMesher<float> aTriangularMesher;

    clock_t start, finish;

    start = clock();

    aTriangularMesher.remesh(anEdgeMesh, meshSize);

    try 
    {
        std::vector<CGeoCoordinate<float>> sInteriorPoints;
        aTriangularMesher.generate(anEdgeMesh, maxEle, sInteriorPoints, meshSize, meshSize, meshSize, tol);
    } 
    catch (std::vector<SMshAdvancingFrontEdge<float>>& advFront)
    {
        CMshMesh<float> aMesh;

		aSurfaceMesh = aTriangularMesher.mesh();

        for (int i = 0; i < aSurfaceMesh.nbNodes(); i++) {
            int aNodeId = aSurfaceMesh.nodeId(i);
            CMshNode<float> aNode = aSurfaceMesh.node(aNodeId);

            aMesh.addNode(aNodeId, aNode);
        }

        for (int i = 0; i < advFront.size(); i++) {

            if (advFront[i].remove)
                continue;

            CMshElement<float> anElement(ET_BEAM);
            anElement.addNodeId(advFront[i].nodeId[0]);
            anElement.addNodeId(advFront[i].nodeId[1]);
            aMesh.addElement(aMesh.nextElementId(), anElement);
        }

	    T.setMesh(aMesh);
        aPosGmsh.save(T, "adv_front.msh", "beams");
    }

    aSurfaceMesh = aTriangularMesher.mesh();

    for (int i = 0; i < 3; i++)
    {
        aTriangularMesher.remesh(aSurfaceMesh, meshSize, tol);
        for (int j = 0; j < 2; j++)
        {
            aTriangularMesher.flipEdges(aSurfaceMesh, tol);
            aTriangularMesher.relaxNodes(aSurfaceMesh, tol);
        }
    }

    finish = clock();
    std::cout << "Finished in about " << std::setprecision(2) << (float)(finish - start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    std::cout << "Number of elements: " << aSurfaceMesh.nbElements() << std::endl;

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "surface.msh", "tris");
}

int main(int argc, char* argv[])
{
    //CMshMesh<float> anEdgeMesh = GenerateBoundary(1E-3);
    CMshMesh<float> anEdgeMesh = LoadBoundary("_edge.msh", 1E-3);

    float m = 1.0;
    for (int i = 0; i < 200; i++)
    {
        std::cout << "i = " << i << std::endl;
        std::cout << "Mesh size = " << std::setprecision(4) << m << std::endl;
        GenerateMesh(anEdgeMesh, m, 999, 0.02);
        m += 0.01;
    }
}
