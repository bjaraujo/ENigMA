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

#include "FemTriangle.hpp"
#include "MshBasicMesher.hpp"
#include "SleSystem.hpp"
#include "PdeEquation.hpp"
#include "PdeBoundaryCondition.hpp"
#include "MatMaterial.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::fem;
using namespace ENigMA::mesh;
using namespace ENigMA::material;
using namespace ENigMA::pde;
using namespace ENigMA::post;

void runExample56()
{

/***************************************************************
    Introduction to finite elements in engineering, 3rd Edition,
    Tirupathi R. Chandrupatla, Ashok D. Belegundu
    Example 5.6, Page 148
****************************************************************/

    std::cout << "*************** Example 5.6 ***************" << std::endl;

    CMshMesh<double> aMesh;

    CMshNode<double> aNode1(3.0, 0.0, 0.0);
    CMshNode<double> aNode2(3.0, 2.0, 0.0);
    CMshNode<double> aNode3(0.0, 2.0, 0.0);
    CMshNode<double> aNode4(0.0, 0.0, 0.0);

    aMesh.addNode(1, aNode1);
    aMesh.addNode(2, aNode2);
    aMesh.addNode(3, aNode3);
    aMesh.addNode(4, aNode4);

    CMshElement<double> anElement1;

    anElement1.setElementType(ET_TRIANGLE);
    anElement1.setThickness(0.5);
    anElement1.addNodeId(1);
    anElement1.addNodeId(2);
    anElement1.addNodeId(4);

    aMesh.addElement(1, anElement1);

    CMshElement<double> anElement2;

    anElement2.setElementType(ET_TRIANGLE);
    anElement2.setThickness(0.5);
    anElement2.addNodeId(3);
    anElement2.addNodeId(4);
    anElement2.addNodeId(2);

    aMesh.addElement(2, anElement2);

    CMatMaterial<double> aMaterial;

    aMaterial.addProperty(PT_ELASTIC_MODULUS, 30E+6);
    aMaterial.addProperty(PT_POISSON_COEFFICIENT, 0.25);

    CPdeField<double> u;

    u.setMesh(aMesh);
    u.setMaterial(aMaterial);
    u.setSimulationType(ST_STRUCTURAL);
    u.setDiscretMethod(DM_FEM);
    u.setDiscretOrder(DO_LINEAR);
    u.setDiscretLocation(DL_NODE);
    u.setNbDofs(2);

    // x
    u.setFixedValue(aMesh.nodeIndex(3), 0, 0.0);
    u.setFixedValue(aMesh.nodeIndex(4), 0, 0.0);

    // y
    u.setFixedValue(aMesh.nodeIndex(1), 1, 0.0);
    u.setFixedValue(aMesh.nodeIndex(3), 1, 0.0);
    u.setFixedValue(aMesh.nodeIndex(4), 1, 0.0);

    u.setSource(aMesh.nodeIndex(2), 1, -1000.0);

    // Stress
    CPdeEquation<double> aPdeEquation(laplacian<double>(u) = 0.0);

    aPdeEquation.setSources(u);

    aPdeEquation.setElimination(u);

    aPdeEquation.solve(u);

    std::cout << u.u << std::endl;

}

void runExample91()
{

/***************************************************************
    Introduction to finite elements in engineering, 3rd Edition,
    Tirupathi R. Chandrupatla, Ashok D. Belegundu
    Example 9.1, Page 284
****************************************************************/

    std::cout << "*************** Example 9.1 ***************" << std::endl;

    CMshMesh<double> aMesh;

    CMshNode<double> aNode1(0.0, 1.0, 1.0);
    CMshNode<double> aNode2(0.0, 0.0, 1.0);
    CMshNode<double> aNode3(1.0, 0.0, 1.0);
    CMshNode<double> aNode4(0.0, 0.0, 0.0);

    aMesh.addNode(1, aNode1);
    aMesh.addNode(2, aNode2);
    aMesh.addNode(3, aNode3);
    aMesh.addNode(4, aNode4);

    CMshElement<double> anElement;

    anElement.setElementType(ET_TETRAHEDRON);
    anElement.addNodeId(1);
    anElement.addNodeId(2);
    anElement.addNodeId(3);
    anElement.addNodeId(4);

    aMesh.addElement(1, anElement);

    CMatMaterial<double> aMaterial;

    aMaterial.addProperty(PT_ELASTIC_MODULUS, 30E+6);
    aMaterial.addProperty(PT_POISSON_COEFFICIENT, 0.3);

    CPdeField<double> u;

    u.setMesh(aMesh);
    u.setMaterial(aMaterial);
    u.setSimulationType(ST_STRUCTURAL);
    u.setDiscretMethod(DM_FEM);
    u.setDiscretOrder(DO_LINEAR);
    u.setDiscretLocation(DL_NODE);
    u.setNbDofs(3);

    // x
    u.setFixedValue(aMesh.nodeIndex(2), 0, 0.0);
    u.setFixedValue(aMesh.nodeIndex(3), 0, 0.0);
    u.setFixedValue(aMesh.nodeIndex(4), 0, 0.0);

    // y
    u.setFixedValue(aMesh.nodeIndex(2), 1, 0.0);
    u.setFixedValue(aMesh.nodeIndex(3), 1, 0.0);
    u.setFixedValue(aMesh.nodeIndex(4), 1, 0.0);

    // z
    u.setFixedValue(aMesh.nodeIndex(2), 2, 0.0);
    u.setFixedValue(aMesh.nodeIndex(3), 2, 0.0);
    u.setFixedValue(aMesh.nodeIndex(4), 2, 0.0);

    u.setSource(aMesh.nodeIndex(1), 2, -1000.0);

    // Stress
    CPdeEquation<double> aPdeEquation(laplacian<double>(u) = 0.0);

    aPdeEquation.setSources(u);

    aPdeEquation.setElimination(u);

    aPdeEquation.solve(u);

    std::cout << u.u << std::endl;

}

void runExample104()
{

/***************************************************************
    Introduction to finite elements in engineering, 3rd Edition,
    Tirupathi R. Chandrupatla, Ashok D. Belegundu
    Example 10.4, Page 326
****************************************************************/

    std::cout << "*************** Example 10.4 ***************" << std::endl;

    CMshMesh<double> aMesh;

    // Add nodes
    CMshNode<double> aNode1(0.0, 0.00, 0.0);
    CMshNode<double> aNode2(0.4, 0.00, 0.0);
    CMshNode<double> aNode3(0.4, 0.15, 0.0);
    CMshNode<double> aNode4(0.4, 0.30, 0.0);
    CMshNode<double> aNode5(0.0, 0.30, 0.0);

    aMesh.addNode(1, aNode1);
    aMesh.addNode(2, aNode2);
    aMesh.addNode(3, aNode3);
    aMesh.addNode(4, aNode4);
    aMesh.addNode(5, aNode5);

    // Add triangle elements
    CMshElement<double> anElement1;

    anElement1.setElementType(ET_TRIANGLE);
    anElement1.setThickness(1.0);
    anElement1.addNodeId(1);
    anElement1.addNodeId(2);
    anElement1.addNodeId(3);

    CMshElement<double> anElement2;

    anElement2.setElementType(ET_TRIANGLE);
    anElement2.setThickness(1.0);
    anElement2.addNodeId(5);
    anElement2.addNodeId(1);
    anElement2.addNodeId(3);

    CMshElement<double> anElement3;

    anElement3.setElementType(ET_TRIANGLE);
    anElement3.setThickness(1.0);
    anElement3.addNodeId(5);
    anElement3.addNodeId(3);
    anElement3.addNodeId(4);

    aMesh.addElement(1, anElement1);
    aMesh.addElement(2, anElement2);
    aMesh.addElement(3, anElement3);

    CMatMaterial<double> aMaterial;

    aMaterial.addProperty(PT_THERMAL_CONDUCTIVITY, 1.5);

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setMaterial(aMaterial);
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    T.setFixedValue(aMesh.nodeIndex(4), 180.0);
    T.setFixedValue(aMesh.nodeIndex(5), 180.0);

    // Heat convection on edge elements
    CPdeBoundaryCondition<double> aConvectiveBC(BT_HEAT_CONVECTIVE);

    aConvectiveBC.addCondition(CT_HEAT_TRANSFER_COEFFICIENT, 50);
    aConvectiveBC.addCondition(CT_HEAT_INFINITESIMAL_TEMPERATURE, 25);
    aConvectiveBC.setLocation(BL_EDGE);

    T.addBCElement(1, 1, aConvectiveBC);    // bc 1, tri 1, edge 1
    T.addBCElement(3, 1, aConvectiveBC);    // bc 2, tri 3, edge 1

    // Steady-state temperature equation
    CPdeEquation<double> aPdeEquation(laplacian<double>(T) = 0.0);

    //aPdeEquation.setElimination(T);
    aPdeEquation.setPenaltyFactor(T, 1000.0);

    aPdeEquation.solve(T);

    std::cout << T.u << std::endl;    

}

void runExample105()
{

/***************************************************************
    Introduction to finite elements in engineering, 3rd Edition,
    Tirupathi R. Chandrupatla, Ashok D. Belegundu
    Example 10.5, Page 334
****************************************************************/

    std::cout << "*************** Example 10.5 ***************" << std::endl;

    CMshMesh<double> aMesh;

    CMshNode<double> aNode1(0.0, 0.0, 0.0);
    CMshNode<double> aNode2(2.0, 1.5, 0.0);
    CMshNode<double> aNode3(4.0, 0.0, 0.0);
    CMshNode<double> aNode4(4.0, 3.0, 0.0);
    CMshNode<double> aNode5(0.0, 3.0, 0.0);

    aMesh.addNode(1, aNode1);
    aMesh.addNode(2, aNode2);
    aMesh.addNode(3, aNode3);
    aMesh.addNode(4, aNode4);
    aMesh.addNode(5, aNode5);

    CMshElement<double> anElement1;

    anElement1.setElementType(ET_TRIANGLE);
    anElement1.addNodeId(1);
    anElement1.addNodeId(3);
    anElement1.addNodeId(2);

    aMesh.addElement(1, anElement1);

    CMshElement<double> anElement2;

    anElement2.setElementType(ET_TRIANGLE);
    anElement2.addNodeId(3);
    anElement2.addNodeId(4);
    anElement2.addNodeId(2);

    aMesh.addElement(2, anElement2);

    CMshElement<double> anElement3;

    anElement3.setElementType(ET_TRIANGLE);
    anElement3.addNodeId(4);
    anElement3.addNodeId(5);
    anElement3.addNodeId(2);

    aMesh.addElement(3, anElement3);

    CMshElement<double> anElement4;

    anElement4.setElementType(ET_TRIANGLE);
    anElement4.addNodeId(5);
    anElement4.addNodeId(1);
    anElement4.addNodeId(2);

    aMesh.addElement(4, anElement4);

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FEM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_GENERIC);
    T.setNbDofs(1);

    T.setFixedValue(aMesh.nodeIndex(3), 0.0);
    T.setFixedValue(aMesh.nodeIndex(4), 0.0);
    T.setFixedValue(aMesh.nodeIndex(5), 0.0);

    // Torsion
    CPdeEquation<double> aPdeEquation(laplacian<double>(T) = source<double>(T, 2.0));

    aPdeEquation.setElimination(T);

    aPdeEquation.solve(T);

    //std::cout << T.u << std::endl;

    double c = 0.0;

    for (Integer i = 0; i < aMesh.nbElements(); ++i)
    {

        Integer anElementId = aMesh.elementId(i);

        CMshElement<double> anElement = aMesh.element(anElementId);

        double sum = 0.0;

        for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
        {

            sum += T.u[aMesh.nodeIndex(anElement.nodeId(j))];

        }

        double area = 3.0;

        c += 8.0 * area / 3.0 * sum;

    }

    std::cout << "alpha = " << 1.0 / c << std::endl;
    
}

int main(int argc, char** argv)
{

    runExample56();
    runExample91();
    runExample104();
    runExample105();

}
