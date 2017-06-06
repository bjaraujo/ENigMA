// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

// QT includes
#include <QApplication>

#include "FemFlow.h"

extern int qInitResources_icons();

int main(int argc, char** argv)
{

    // QT Stuff
    QApplication app(argc, argv);

    qInitResources_icons();

    FemFlow myFemFlow;
    myFemFlow.show();

    myFemFlow.init();

    return app.exec();

}

/*
#include <iostream>
#include <iomanip>

#include "FemTriangle.hpp"
#include "FemFlowTriangle.hpp"
#include "MshBasicMesher.hpp"
#include "MshTriangleMesher.hpp"
#include "SleSystem.hpp"
#include "PdeEquation.hpp"
#include "PosGnuplot.hpp"
#include "PosGmsh.hpp"
#include "AnaTemperature.hpp"

using namespace Eigen;

using namespace ENigMA::fem;
using namespace ENigMA::fem::flow;
using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::analytical;

void flow2D_tri()
{

std::cout << "**** Two-dimensional flow (tri) ****" << std::endl;

CGeoCoordinate<double> aVertex1(0.0, 0.0, 0.0);
CGeoCoordinate<double> aVertex2(1.0, 0.0, 0.0);
CGeoCoordinate<double> aVertex3(1.0, 1.0, 0.0);
CGeoCoordinate<double> aVertex4(0.0, 1.0, 0.0);

CGeoQuadrilateral<double> aQuadrilateral;

aQuadrilateral.addVertex(aVertex1);
aQuadrilateral.addVertex(aVertex2);
aQuadrilateral.addVertex(aVertex3);
aQuadrilateral.addVertex(aVertex4);

CMshBasicMesher<double> aBasicMesher;

CMshMesh<double> aMesh = aBasicMesher.meshGeometry(aQuadrilateral, 25, 25, true);

for (Integer i = 0; i < aMesh.nbElements(); ++i)
{

Integer anElementId = aMesh.elementId(i);
aMesh.element(anElementId).setThickness(1.0);

}

CPdeField<double> u, v, p;

u.setMesh(aMesh);
u.setDiscretMethod(DM_FEM);
u.setDiscretOrder(DO_LINEAR);
u.setDiscretLocation(DL_NODE);
u.setSimulationType(ST_FLOW);
u.setNbDofs(1);

v.setMesh(aMesh);
v.setDiscretMethod(DM_FEM);
v.setDiscretOrder(DO_LINEAR);
v.setDiscretLocation(DL_NODE);
v.setSimulationType(ST_FLOW);
v.setNbDofs(1);

p.setMesh(aMesh);
p.setDiscretMethod(DM_FEM);
p.setDiscretOrder(DO_LINEAR);
p.setDiscretLocation(DL_NODE);
p.setSimulationType(ST_FLOW);
p.setNbDofs(1);

// Set BC and initial conditions

u.u.resize(u.mesh().nbNodes());
v.u.resize(v.mesh().nbNodes());
p.u.resize(p.mesh().nbNodes());

double Ulid = 1.0;

for (Integer i = 0; i < u.mesh().nbNodes(); ++i)
{

Integer aNodeId = u.mesh().nodeId(i);
CMshNode<double> aNode = u.mesh().node(aNodeId);

if (fabs(aNode.x() - 0.0) < 1E-6)
u.setFixedValue(i, 0.0);

if (fabs(aNode.x() - 1.0) < 1E-6)
u.setFixedValue(i, 0.0);

if (fabs(aNode.y() - 0.0) < 1E-6)
u.setFixedValue(i, 0.0);

if (fabs(aNode.y() - 1.0) < 1E-6)
u.setFixedValue(i, Ulid);

u.u(i) = 0.0;

}

for (Integer i = 0; i < v.mesh().nbNodes(); ++i)
{

Integer aNodeId = v.mesh().nodeId(i);
CMshNode<double> aNode = v.mesh().node(aNodeId);

if (fabs(aNode.x() - 0.0) < 1E-6)
v.setFixedValue(i, 0.0);

if (fabs(aNode.x() - 1.0) < 1E-6)
v.setFixedValue(i, 0.0);

if (fabs(aNode.y() - 0.0) < 1E-6)
v.setFixedValue(i, 0.0);

if (fabs(aNode.y() - 1.0) < 1E-6)
v.setFixedValue(i, 0.0);

v.u(i) = 0.0;

}

for (Integer i = 0; i < p.mesh().nbNodes(); ++i)
{

Integer aNodeId = p.mesh().nodeId(i);
CMshNode<double> aNode = p.mesh().node(aNodeId);

if (fabs(aNode.x() - 0.0) < 1E-6 && fabs(aNode.y() - 0.0) < 1E-6)
p.setFixedValue(i, 0.0);

p.u(i) = 0.0;

}

double mu = 1.0;  // dynamic viscosity
double rho = 1.0; // density

double nu = mu / rho; // kinematic viscosity

std::cout << "Re = " << Ulid / nu << std::endl;

SparseMatrix<double> Mu = ddt<double>(u).matrixA;
SparseMatrix<double> Mv = ddt<double>(v).matrixA;

SparseMatrix<double> G1 = gradient<double>(u, CP_X).matrixA;
SparseMatrix<double> G2 = gradient<double>(v, CP_Y).matrixA;

double dt = 1E-4;
Integer nIter = 10;

// Flow in a rectangle
for (Integer i = 0; i < nIter; ++i)
{

std::cout << "Iteration = " << (i + 1) << std::endl;
std::cout << "Time = " << dt * (i + 1) << std::endl;

SparseMatrix<double> Ku = laplacian<double>(u).matrixA;
SparseMatrix<double> Kv = laplacian<double>(v).matrixA;

SparseMatrix<double> C = divergence<double>(u, v, dt).matrixA;

// Velocity calculation
CPdeEquation<double> aPdeEquationU(ddt<double>(u) = dt * (nu * Ku - C) * u.u);
aPdeEquationU.setElimination(u);
aPdeEquationU.solve(u);

CPdeEquation<double> aPdeEquationV(ddt<double>(v) = dt * (nu * Kv - C) * v.u);
aPdeEquationV.setElimination(v);
aPdeEquationV.solve(v);

// Pressure calculation
CPdeEquation<double> aPdeEquationP(laplacian<double>(p) = rho / dt * (G1 * u.u + G2 * v.u));
aPdeEquationP.setElimination(p);
aPdeEquationP.solve(p);

// Velocity correction
CPdeEquation<double> aPdeEquationUCor(ddt<double>(u) = -dt / rho * G1 * p.u);
CPdeEquation<double> aPdeEquationVCor(ddt<double>(v) = -dt / rho * G2 * p.u);

aPdeEquationUCor.setElimination(u);
aPdeEquationVCor.setElimination(v);

aPdeEquationUCor.solve(u);
aPdeEquationVCor.solve(v);

}

CPosGmsh<double> aPosGmsh;

CPdeField<double> U;

U.setMesh(aMesh);
U.setDiscretLocation(DL_NODE);
U.setNbDofs(3);

U.u.resize(aMesh.nbNodes() * 3, 1);

for (Integer i = 0; i < aMesh.nbNodes(); ++i)
{

Integer aNodeId = aMesh.nodeId(i);

U.u(i * 3 + 0) = u.u[aNodeId];
U.u(i * 3 + 1) = v.u[aNodeId];
U.u(i * 3 + 2) = 0.0;

}

aPosGmsh.save(U, "fem_flow_tri_u_2D.msh", "Velocity u");
aPosGmsh.save(p, "fem_flow_tri_p_2D.msh", "Pressure p");

}

void flow2D_tetra()
{

std::cout << "**** Two-dimensional flow (tetra) ****" << std::endl;

CGeoCoordinate<double> aVertex1(0.0, 0.0, -1.0);
CGeoCoordinate<double> aVertex2(1.0, 0.0, -1.0);
CGeoCoordinate<double> aVertex3(1.0, 1.0, -1.0);
CGeoCoordinate<double> aVertex4(0.0, 1.0, -1.0);
CGeoCoordinate<double> aVertex5(0.0, 0.0, +1.0);
CGeoCoordinate<double> aVertex6(1.0, 0.0, +1.0);
CGeoCoordinate<double> aVertex7(1.0, 1.0, +1.0);
CGeoCoordinate<double> aVertex8(0.0, 1.0, +1.0);

CGeoHexahedron<double> aHexahedron;

aHexahedron.addVertex(aVertex1);
aHexahedron.addVertex(aVertex2);
aHexahedron.addVertex(aVertex3);
aHexahedron.addVertex(aVertex4);
aHexahedron.addVertex(aVertex5);
aHexahedron.addVertex(aVertex6);
aHexahedron.addVertex(aVertex7);
aHexahedron.addVertex(aVertex8);

CMshBasicMesher<double> aBasicMesher;

CMshMesh<double> aMesh = aBasicMesher.meshGeometry(aHexahedron, 25, 25, 1, true);

CPdeField<double> u, v, w, p;

u.setMesh(aMesh);
u.setDiscretMethod(DM_FEM);
u.setDiscretOrder(DO_LINEAR);
u.setDiscretLocation(DL_NODE);
u.setSimulationType(ST_FLOW);
u.setNbDofs(1);

v.setMesh(aMesh);
v.setDiscretMethod(DM_FEM);
v.setDiscretOrder(DO_LINEAR);
v.setDiscretLocation(DL_NODE);
v.setSimulationType(ST_FLOW);
v.setNbDofs(1);

w.setMesh(aMesh);
w.setDiscretMethod(DM_FEM);
w.setDiscretOrder(DO_LINEAR);
w.setDiscretLocation(DL_NODE);
w.setSimulationType(ST_FLOW);
w.setNbDofs(1);

p.setMesh(aMesh);
p.setDiscretMethod(DM_FEM);
p.setDiscretOrder(DO_LINEAR);
p.setDiscretLocation(DL_NODE);
p.setSimulationType(ST_FLOW);
p.setNbDofs(1);

// Set BC and initial conditions

u.u.resize(u.mesh().nbNodes());
v.u.resize(v.mesh().nbNodes());
w.u.resize(w.mesh().nbNodes());
p.u.resize(p.mesh().nbNodes());

double Ulid = 1.0;

for (Integer i = 0; i < u.mesh().nbNodes(); ++i)
{

Integer aNodeId = u.mesh().nodeId(i);
CMshNode<double> aNode = u.mesh().node(aNodeId);

if (fabs(aNode.x() - 0.0) < 1E-6)
u.setFixedValue(i, 0.0);

if (fabs(aNode.x() - 1.0) < 1E-6)
u.setFixedValue(i, 0.0);

if (fabs(aNode.y() - 0.0) < 1E-6)
u.setFixedValue(i, 0.0);

if (fabs(aNode.y() - 1.0) < 1E-6)
u.setFixedValue(i, Ulid);

u.u(i) = 0.0;

}

for (Integer i = 0; i < v.mesh().nbNodes(); ++i)
{

Integer aNodeId = v.mesh().nodeId(i);
CMshNode<double> aNode = v.mesh().node(aNodeId);

if (fabs(aNode.x() - 0.0) < 1E-6)
v.setFixedValue(i, 0.0);

if (fabs(aNode.x() - 1.0) < 1E-6)
v.setFixedValue(i, 0.0);

if (fabs(aNode.y() - 0.0) < 1E-6)
v.setFixedValue(i, 0.0);

if (fabs(aNode.y() - 1.0) < 1E-6)
v.setFixedValue(i, 0.0);

v.u(i) = 0.0;

}

for (Integer i = 0; i < w.mesh().nbNodes(); ++i)
{

Integer aNodeId = w.mesh().nodeId(i);
CMshNode<double> aNode = w.mesh().node(aNodeId);

w.u(i) = 0.0;

}

for (Integer i = 0; i < p.mesh().nbNodes(); ++i)
{

Integer aNodeId = p.mesh().nodeId(i);
CMshNode<double> aNode = p.mesh().node(aNodeId);

if (fabs(aNode.x() - 0.0) < 1E-6 && fabs(aNode.y() - 0.0) < 1E-6)
p.setFixedValue(i, 0.0);

p.u(i) = 0.0;

}

double mu = 1.0;  // dynamic viscosity
double rho = 1.0; // density

double nu = mu / rho; // kinematic viscosity

std::cout << "Re = " << Ulid / nu << std::endl;

SparseMatrix<double> G1 = gradient<double>(u, CP_X).matrixA;
SparseMatrix<double> G2 = gradient<double>(v, CP_Y).matrixA;

double dt = 1E-4;
Integer nIter = 10;

// Flow in a box
for (Integer i = 0; i < nIter; ++i)
{

std::cout << "Iteration = " << (i + 1) << std::endl;
std::cout << "Time = " << dt * (i + 1) << std::endl;

SparseMatrix<double> Ku = laplacian<double>(u).matrixA;
SparseMatrix<double> Kv = laplacian<double>(v).matrixA;

SparseMatrix<double> C = divergence<double>(u, v, w, dt).matrixA;

// Velocity calculation
CPdeEquation<double> aPdeEquationU(ddt<double>(u) = dt * (nu * Ku - C) * u.u);
aPdeEquationU.setElimination(u);
aPdeEquationU.solve(u);

CPdeEquation<double> aPdeEquationV(ddt<double>(v) = dt * (nu * Kv - C) * v.u);
aPdeEquationV.setElimination(v);
aPdeEquationV.solve(v);

// Pressure calculation
CPdeEquation<double> aPdeEquationP(laplacian<double>(p) = rho / dt * (G1 * u.u + G2 * v.u));
aPdeEquationP.setElimination(p);
aPdeEquationP.solve(p);

// Velocity correction
CPdeEquation<double> aPdeEquationUCor(ddt<double>(u) = -dt / rho * G1 * p.u);
CPdeEquation<double> aPdeEquationVCor(ddt<double>(v) = -dt / rho * G2 * p.u);

aPdeEquationUCor.setElimination(u);
aPdeEquationVCor.setElimination(v);

aPdeEquationUCor.solve(u);
aPdeEquationVCor.solve(v);

}

CPosGmsh<double> aPosGmsh;

CPdeField<double> U;

U.setMesh(aMesh);
U.setDiscretLocation(DL_NODE);
U.setNbDofs(3);

U.u.resize(aMesh.nbNodes() * 3, 1);

for (Integer i = 0; i < aMesh.nbNodes(); ++i)
{

Integer aNodeId = aMesh.nodeId(i);

U.u(i * 3 + 0) = u.u[aNodeId];
U.u(i * 3 + 1) = v.u[aNodeId];
U.u(i * 3 + 2) = w.u[aNodeId];

}

aPosGmsh.save(U, "fem_flow_tetra_u_2D.msh", "Velocity U");
aPosGmsh.save(p, "fem_flow_tetra_p_2D.msh", "Pressure p");

}

int main(int argc, char** argv)
{

flow2D_tri();
//flow2D_tetra();

}
*/
