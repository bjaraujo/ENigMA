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

void GenerateMesh(const double meshSize, const int maxEle, const double tol)
{
    CPdeField<double> T;
    CPosGmsh<double> aPosGmsh;
    CMshMesh<double> aSurfaceMesh;

    //CMshQuadrilateralMesher<double> aQuadrilateralMesher;
    CMshTriangleMesher<double> aTriangularMesher;

    aPosGmsh.load(T, "_edge947.msh");
    CMshMesh<double> anEdgeMesh = T.mesh();

    anEdgeMesh.mergeNodes(1E-2);
    anEdgeMesh.generateFaces(1E-2);

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "edge_quads.msh", "beams");

    clock_t start, finish;

    start = clock();

    //aTriangularMesher.remesh(anEdgeMesh, meshSize);

    try {
        aTriangularMesher.generate(anEdgeMesh, maxEle, meshSize, 0.1 * meshSize, 10.0 * meshSize, tol);
    } catch (std::vector<SMshTriangleAdvancingFrontEdge<double>>& advFront) {

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

    finish = clock();

    std::cout << "Finished in about " << std::setprecision(2) << (double)(finish - start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    aSurfaceMesh = aTriangularMesher.mesh();

    std::cout << "Number of elements: " << aSurfaceMesh.nbElements() << std::endl;

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "surface_quads.msh", "quads");
}

int main(int argc, char* argv[])
{
    GenerateMesh(1.0, 999, 1E-2);
}
