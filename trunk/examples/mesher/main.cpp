// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iostream>
#include <iomanip>

#include <time.h>

#include "MshTriangleMesher.hpp"
#include "MshQuadrilateralMesher.hpp"
#include "PosGmsh.hpp"
#include "StlUtils.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::pde;
using namespace ENigMA::stl;

void Example1(const double meshSize, const int maxEle, const double tol)
{

    CMshTriangleMesher<double> aTriangleMesher;

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

    anEdgeMesh.generateFaces(1E-3);

    clock_t start, finish;

    start = clock();

    aTriangleMesher.remesh(anEdgeMesh, meshSize);
    aTriangleMesher.generate(anEdgeMesh, maxEle, meshSize, 0.1, tol);

    aTriangleMesher.collapseEdges(meshSize, tol);
    
    aTriangleMesher.relaxNodes(tol);

    finish = clock();

    std::cout << "Finished in about " << std::setprecision(2) << (double)(finish - start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    CMshMesh<double> aSurfaceMesh;
    aSurfaceMesh = aTriangleMesher.mesh();

    std::cout << "Number of elements: " << aSurfaceMesh.nbElements() << std::endl;

    CPdeField<double> T;
    CPosGmsh<double> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "edge_tris.msh", "beams");

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "surface_tris.msh", "tris");

}

void Example2(const double meshSize, const int maxEle, const double tol)
{

    CMshQuadrilateralMesher<double> aQuadrilateralMesher;

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

    CMshNode<double> aNode5;
    aNode5 << 0.6, 0.5, 0.0;
    anEdgeMesh.addNode(4, aNode5);

    CMshNode<double> aNode6;
    aNode6 << 0.4, 0.5, 0.0;
    anEdgeMesh.addNode(5, aNode6);

    CMshNode<double> aNode7;
    aNode7 << 0.5, 0.7, 0.0;
    anEdgeMesh.addNode(6, aNode7);

    CMshElement<double> anElement5(ET_BEAM);
    anElement5.addNodeId(4);
    anElement5.addNodeId(5);
    anEdgeMesh.addElement(4, anElement5);

    CMshElement<double> anElement6(ET_BEAM);
    anElement6.addNodeId(5);
    anElement6.addNodeId(6);
    anEdgeMesh.addElement(5, anElement6);

    CMshElement<double> anElement7(ET_BEAM);
    anElement7.addNodeId(6);
    anElement7.addNodeId(4);
    anEdgeMesh.addElement(6, anElement7);

    anEdgeMesh.generateFaces(1E-3);

    clock_t start, finish;

    start = clock();

    aQuadrilateralMesher.remesh(anEdgeMesh, meshSize);
    aQuadrilateralMesher.generate(anEdgeMesh, maxEle, meshSize, 0.1, tol);

    finish = clock();

    std::cout << "Finished in about " << std::setprecision(2) << (double)(finish - start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    CMshMesh<double> aSurfaceMesh;
    aSurfaceMesh = aQuadrilateralMesher.mesh();

    std::cout << "Number of elements: " << aSurfaceMesh.nbElements() << std::endl;

    CPdeField<double> T;
    CPosGmsh<double> aPosGmsh;

    T.setMesh(anEdgeMesh);
    aPosGmsh.save(T, "edge_quads.msh", "beams");

    T.setMesh(aSurfaceMesh);
    aPosGmsh.save(T, "surface_quads.msh", "quads");

}

int main(int argc, char *argv[])
{

    if (argc > 1)
    {
        //Example1(atof(argv[1]), atoi(argv[2]), atof(argv[3]));
        Example2(atof(argv[1]), atoi(argv[2]), atof(argv[3]));
    }

}
