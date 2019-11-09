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

#include "BemTriangle.hpp"
#include "MshBasicMesher.hpp"
#include "SleSystem.hpp"
#include "PdeEquation.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::bem;
using namespace ENigMA::mesh;
using namespace ENigMA::post;

void solve1()
{

    std::cout << "Meshing..." << std::endl;

    CGeoCoordinate<double> aVertex01(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex02(1.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex03(1.0, 0.5, 0.0);
    CGeoCoordinate<double> aVertex04(0.0, 0.5, 0.0);
    CGeoCoordinate<double> aVertex05(0.0, 0.0, 0.5);
    CGeoCoordinate<double> aVertex06(1.0, 0.0, 0.5);
    CGeoCoordinate<double> aVertex07(1.0, 0.5, 0.5);
    CGeoCoordinate<double> aVertex08(0.0, 0.5, 0.5);

    CGeoHexahedron<double> aHexahedron1;

    aHexahedron1.addVertex(aVertex01);
    aHexahedron1.addVertex(aVertex02);
    aHexahedron1.addVertex(aVertex03);
    aHexahedron1.addVertex(aVertex04);
    aHexahedron1.addVertex(aVertex05);
    aHexahedron1.addVertex(aVertex06);
    aHexahedron1.addVertex(aVertex07);
    aHexahedron1.addVertex(aVertex08);

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aHexahedron1, 12, 6, 6, true);

    CMshMesh<double> aMesh1 = aBasicMesher.mesh();

    CGeoCoordinate<double> aVertex09(0.1, 0.1, 0.39);
    CGeoCoordinate<double> aVertex10(0.9, 0.1, 0.39);
    CGeoCoordinate<double> aVertex11(0.9, 0.4, 0.39);
    CGeoCoordinate<double> aVertex12(0.1, 0.4, 0.39);
    CGeoCoordinate<double> aVertex13(0.1, 0.1, 0.4);
    CGeoCoordinate<double> aVertex14(0.9, 0.1, 0.4);
    CGeoCoordinate<double> aVertex15(0.9, 0.4, 0.4);
    CGeoCoordinate<double> aVertex16(0.1, 0.4, 0.4);

    CGeoHexahedron<double> aHexahedron2;

    aHexahedron2.addVertex(aVertex09);
    aHexahedron2.addVertex(aVertex10);
    aHexahedron2.addVertex(aVertex11);
    aHexahedron2.addVertex(aVertex12);
    aHexahedron2.addVertex(aVertex13);
    aHexahedron2.addVertex(aVertex14);
    aHexahedron2.addVertex(aVertex15);
    aHexahedron2.addVertex(aVertex16);

    aBasicMesher.generate(aHexahedron2, 11, 5, 1, true);

    CMshMesh<double> aMesh2 = aBasicMesher.mesh();

    std::cout << "Extracting surface mesh..." << std::endl;

    CMshMesh<double> aSurfaceMesh1 = aMesh1.extractBoundary(1E-12);
    CMshMesh<double> aSurfaceMesh2 = aMesh2.extractBoundary(1E-12);

    aSurfaceMesh1.invert();

    CMshMesh<double> aSurfaceMesh;

    aSurfaceMesh.addMesh(aSurfaceMesh1);
    aSurfaceMesh.addMesh(aSurfaceMesh2);

    std::vector<CBemTriangle<double> > aBemMesh;

    for (Integer i = 0; i < aSurfaceMesh.nbElements(); ++i)
    {

        Integer anElementId = aSurfaceMesh.elementId(i);

        CMshElement<double> anElement = aSurfaceMesh.element(anElementId);

        if (anElement.elementType() != ET_TRIANGLE)
            continue;

        CBemTriangle<double> aTriangle;

        for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
        {

            Integer aNodeId = anElement.nodeId(i);

            CMshNode<double> aNode = aSurfaceMesh.node(aNodeId);

            CGeoCoordinate<double> aVertex = aNode;

            aTriangle.addVertex(aVertex);

        }

        aBemMesh.push_back(aTriangle);

    }

    // Set boundary conditions

    std::cout << "Setting boundary conditions..." << std::endl;

    std::vector<double> q;
    std::vector<double> Tf;
    std::vector<bool> Tfixed;

    for (Integer i = 0; i < static_cast<Integer>(aBemMesh.size()); ++i)
    {

        Tf.push_back(0.0);
        Tfixed.push_back(false);
        q.push_back(0.0);

        aBemMesh[i].calculateCentroid();

        if (fabs(aBemMesh[i].centroid().x() - 0.0) < 1E-12)
        {
            Tf[i] = 0.0;
            Tfixed[i] = true;
        }

        if (fabs(aBemMesh[i].centroid().x() - 1.0) < 1E-12)
        {
            Tf[i] = 1.0;
            Tfixed[i] = true;
        }

    }

    // Solve
    std::cout << "Assembling system..." << std::endl;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A;
    Eigen::Matrix<double, Eigen::Dynamic, 1> b;
    Eigen::Matrix<double, Eigen::Dynamic, 1> T;

    std::cout << aBemMesh.size() << " x " << aBemMesh.size() << std::endl;

    A.resize(aBemMesh.size(), aBemMesh.size());
    A.setZero();

    b.resize(aBemMesh.size());
    b.setZero();

    const double hcoeff = 1.0;
    const double thcond = 1.0;

    double perc;

    for (Integer i = 0; i < static_cast<Integer>(aBemMesh.size()); ++i)
    {

        perc = (double) i / aBemMesh.size() * 100;

        if (i % 50 == 0)
        {
            std::cout << perc << " %..." << std::endl;
        }

        for (Integer j = 0; j < static_cast<Integer>(aBemMesh.size()); ++j)
        {

            double hij, gij;

            aBemMesh[i].laplacianCoeff(i, j, aBemMesh[j], hij, gij);

            double Delta = aBemMesh[j].area();

            //std::cout << "Delta = " << Delta << std::endl;

            /*
            if (i >= aBemMesh.size() - 144*2)
            {
                hij *= 0.5;
                gij *= 0.5;
                q[j] = 1;
            }
            */

            // H
            if (i == j)
                A(i, j) += hij - 0.5;
            else
                A(i, j) += hij;

            if (Tfixed[j])
            {
                // H^
                A(i, j) += -hcoeff / thcond / Delta * gij;

                b[i] += -hcoeff / thcond / Delta * gij * Tf[j];
            }
            else
            {
                b[i] += -q[j] / thcond * gij;
            }

        }

    }
    
    std::cout << "Solving system..." << std::endl;

    T = A.lu().solve(b);

    CPdeField<double> Ti;
    CPosGmsh<double> aPosGmsh;

    Ti.setMesh(aSurfaceMesh);
    Ti.setDiscretLocation(DL_ELEMENT_CENTER);
    Ti.setNbDofs(1);

    Ti.u = T;

    std::cout << "Saving results..." << std::endl;

    aPosGmsh.save(Ti, "bem_heat.pos", "heat");

    std::cout << "Done." << std::endl;

}

void solve2()
{

    CGeoCoordinate<double> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex3(1.0, 0.5, 0.0);
    CGeoCoordinate<double> aVertex4(0.0, 0.5, 0.0);
    CGeoCoordinate<double> aVertex5(0.0, 0.0, 0.5);
    CGeoCoordinate<double> aVertex6(1.0, 0.0, 0.5);
    CGeoCoordinate<double> aVertex7(1.0, 0.5, 0.5);
    CGeoCoordinate<double> aVertex8(0.0, 0.5, 0.5);

    CGeoHexahedron<double> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    const Integer nu = 12;
    const Integer nv = 6;
    const Integer nw = 6;

    std::cout << "Meshing..." << std::endl;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    std::cout << "Extracting surface mesh..." << std::endl;

    CMshMesh<double> aSurfaceMesh = aMesh.extractBoundary(1E-12);

    // Invert elements
    for (Integer i = 0; i < aSurfaceMesh.nbElements(); ++i)
    {

        Integer anElementId = aSurfaceMesh.elementId(i);
        
        aSurfaceMesh.element(anElementId).invert();

    }

    CMshNode<double> aNode1(0.5, 0.25, 0.25);
    CMshNode<double> aNode2(0.6, 0.25, 0.25);
    CMshNode<double> aNode3(0.6, 0.35, 0.25);

    Integer aNodeId = aSurfaceMesh.nextNodeId();

    aSurfaceMesh.addNode(aNodeId + 0, aNode1);
    aSurfaceMesh.addNode(aNodeId + 1, aNode2);
    aSurfaceMesh.addNode(aNodeId + 2, aNode3);

    CMshElement<double> anElement;

    anElement.setElementType(ET_TRIANGLE);

    anElement.addNodeId(aNodeId + 0);
    anElement.addNodeId(aNodeId + 1);
    anElement.addNodeId(aNodeId + 2);

    aSurfaceMesh.addElement(aSurfaceMesh.nextElementId(), anElement);

    anElement.invert();

    aSurfaceMesh.addElement(aSurfaceMesh.nextElementId(), anElement);

    CPdeField<double> Tn;

    Tn.setMesh(aSurfaceMesh);
    Tn.setDiscretMethod(DM_BEM);
    Tn.setDiscretLocation(DL_ELEMENT_CENTER);
    Tn.setDiscretOrder(DO_LINEAR);
    Tn.setSimulationType(ST_THERMAL);
    Tn.setNbDofs(1);

    Tn.mesh().calculateElementCentroid();

    for (Integer i = 0; i < Tn.mesh().nbElements(); ++i)
    {

        Integer anElementId = Tn.mesh().elementId(i);

        CGeoCoordinate<double> aCentroid = Tn.mesh().elementCentroid(anElementId);

        if (fabs(aCentroid.x() - 0.0) < 1E-12)
            Tn.setFixedValue(i,  0.0);

        if (fabs(aCentroid.x() - 1.0) < 1E-12)
            Tn.setFixedValue(i,  1.0);

    }

    std::cout << "Assembling matrix..." << std::endl;

    CPdeEquation<double> aPdeEquation(laplacian<double>(Tn) = 0);

    std::cout << "Solving..." << std::endl;
    aPdeEquation.solve(Tn);

    std::cout << "Saving results..." << std::endl;

    CPosGmsh<double> aPosGmsh;
    aPosGmsh.save(Tn, "bem_pde_heat.pos", "heat");

    std::cout << "Done." << std::endl;

}

int main(int argc, char** argv)
{

    solve1();
    // solve2();

}

