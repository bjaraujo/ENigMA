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

CMshMesh<double> GenerateBoundary(const double tol)
{
    CMshMesh<double> anEdgeMesh;

    CMshNode<double> aNode1;
    aNode1 << 0.0, 0.0, 0.0;
    anEdgeMesh.addNode(0, aNode1);

    CMshNode<double> aNode2;
    aNode2 << 1.0, 0.0, 0.0;
    anEdgeMesh.addNode(1, aNode2);

    CMshNode<double> aNode3;
    aNode3 << 1.0, 1.0, 0.0;
    anEdgeMesh.addNode(2, aNode3);

    CMshNode<double> aNode4;
    aNode4 << 0.0, 1.0, 0.0;
    anEdgeMesh.addNode(3, aNode4);

    CMshElement<double> anElement1(ET_BEAM);
    anElement1.addNodeId(0);
    anElement1.addNodeId(1);
    anEdgeMesh.addElement(0, anElement1);

    CMshElement<double> anElement2(ET_BEAM);
    anElement2.addNodeId(1);
    anElement2.addNodeId(2);
    anEdgeMesh.addElement(1, anElement2);

    CMshElement<double> anElement3(ET_BEAM);
    anElement3.addNodeId(2);
    anElement3.addNodeId(3);
    anEdgeMesh.addElement(2, anElement3);

    CMshElement<double> anElement4(ET_BEAM);
    anElement4.addNodeId(3);
    anElement4.addNodeId(0);
    anEdgeMesh.addElement(3, anElement4);

    anEdgeMesh.generateFaces(tol);
    anEdgeMesh.mergeNodes(tol);

    return anEdgeMesh;
}

CMshMesh<double> LoadBoundary(std::string fileName, const double tol)
{
    CMshMesh<double> anEdgeMesh;

    CPdeField<double> T;
    CPosGmsh<double> aPosGmsh;

    aPosGmsh.load(T, fileName);
    anEdgeMesh = T.mesh();

    anEdgeMesh.mergeNodes(tol);
    anEdgeMesh.generateFaces(tol);

    return anEdgeMesh;
}

void GenerateMesh(CMshMesh<double> anEdgeMesh, const double meshSize, const int maxEle, const double tol)
{
    CMshMesh<double> aSurfaceMesh;

    CPdeField<double> T;
    CPosGmsh<double> aPosGmsh;

    CMshTriangleMesher<double> aTriangularMesher;

    clock_t start, finish;

    start = clock();

    aTriangularMesher.remesh(anEdgeMesh, meshSize);

    try 
    {
        std::vector<CGeoCoordinate<double>> sInteriorPoints;
        aTriangularMesher.generate(anEdgeMesh, maxEle, sInteriorPoints, meshSize, meshSize, tol);
    } 
    catch (std::vector<SMshAdvancingFrontEdge<double>>& advFront) 
    {
        CMshMesh<double> aMesh;

		aSurfaceMesh = aTriangularMesher.mesh();

        for (int i = 0; i < aSurfaceMesh.nbNodes(); i++) {
            int aNodeId = aSurfaceMesh.nodeId(i);
			CMshNode<double> aNode = aSurfaceMesh.node(aNodeId);

            aMesh.addNode(aNodeId, aNode);
        }

        for (int i = 0; i < advFront.size(); i++) {

            if (advFront[i].remove)
                continue;

            CMshElement<double> anElement(ET_BEAM);
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
    std::cout << "Finished in about " << std::setprecision(2) << (double)(finish - start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    std::cout << "Number of elements: " << aSurfaceMesh.nbElements() << std::endl;

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "surface.msh", "tris");
}

int main(int argc, char* argv[])
{
    CMshMesh<double> anEdgeMesh = GenerateBoundary(1E-3);
    //CMshMesh<double> anEdgeMesh = LoadBoundary("_edge.msh", 1E-3);

    GenerateMesh(anEdgeMesh, 0.1, 999, 1E-3);
}
