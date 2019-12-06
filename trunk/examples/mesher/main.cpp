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

void addCircle(CMshMesh<double>& anEdgeMesh, double cx, double cy, double radius)
{
    int aFirstNodeId = anEdgeMesh.nextNodeId();

    for (int i = 0; i < 16; i++) {
        double teta = i / 16.0 * 2 * 3.14159265359;

        CMshNode<double> aNode5(cx - radius * cos(teta), cy + radius * sin(teta), 0.0);
        anEdgeMesh.addNode(anEdgeMesh.nextNodeId(), aNode5);

        if (i > 0) {
            CMshElement<double> anElement(ET_BEAM);
            anElement.addNodeId(anEdgeMesh.nextNodeId() - 2);
            anElement.addNodeId(anEdgeMesh.nextNodeId() - 1);
            anEdgeMesh.addElement(anEdgeMesh.nextElementId(), anElement);
        }
    }

    {
        CMshElement<double> anElement(ET_BEAM);
        anElement.addNodeId(anEdgeMesh.nextNodeId() - 1);
        anElement.addNodeId(aFirstNodeId);
        anEdgeMesh.addElement(anEdgeMesh.nextElementId(), anElement);
    }
}

void GenerateMesh(const double meshSize, const int maxEle, const double tol)
{

    CPdeField<double> T;
    CPosGmsh<double> aPosGmsh;
    CMshMesh<double> aSurfaceMesh;

    CMshQuadrilateralMesher<double> aQuadrilateralMesher;

    CMshMesh<double> anEdgeMesh;

    CMshNode<double> aNode1;
    aNode1 << 0.0, 0.0, 0.0;
    anEdgeMesh.addNode(0, aNode1);

    CMshNode<double> aNode2;
    aNode2 << 400.0, 0.0, 0.0;
    anEdgeMesh.addNode(1, aNode2);

    CMshNode<double> aNode3;
    aNode3 << 400.0, 400.0, 0.0;
    anEdgeMesh.addNode(2, aNode3);

    CMshNode<double> aNode4;
    aNode4 << 0.0, 400.0, 0.0;
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

    addCircle(anEdgeMesh, 200, 200, 100);

    anEdgeMesh.generateFaces(1E-3);

    clock_t start, finish;

    start = clock();

    aQuadrilateralMesher.remesh(anEdgeMesh, meshSize);

    try {
        aQuadrilateralMesher.generate(anEdgeMesh, maxEle, meshSize, 0.1, tol);
    } catch (std::vector<SMshQuadrilateralAdvancingFrontEdge<double>> advFront) {

        CMshMesh<double> aMesh;

		aSurfaceMesh = aQuadrilateralMesher.mesh();

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

    finish = clock();

    std::cout << "Finished in about " << std::setprecision(2) << (double)(finish - start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    aSurfaceMesh = aQuadrilateralMesher.mesh();

    std::cout << "Number of elements: " << aSurfaceMesh.nbElements() << std::endl;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "edge_quads.msh", "beams");

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "surface_quads.msh", "quads");
}

int main(int argc, char* argv[])
{
    GenerateMesh(20, 499, 1E-4);
}
