// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "ui_FemFlow.h"
#include "FemFlow.h"

#include <vtkEventQtSlotConnect.h>
#include <vtkDataObjectToTable.h>
#include <vtkQtTableView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkSTLReader.h>
#include <vtkFeatureEdges.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkAbstractPicker.h>
#include <vtkBandedPolyDataContourFilter.h>
#include <vtkClipPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkContourFilter.h>
#include <vtkCutter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkMarchingCubes.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkConeSource.h>
#include <vtkCellPicker.h>
#include <vtkPlane.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkCamera.h>

#include "FemTriangle.hpp"
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

// Constructor
FemFlow::FemFlow()
{

    this->ui = new Ui_FemFlow;
    this->ui->setupUi(this);

    this->setWindowTitle("Fem Flow");

    // VTK Renderer
    this->m_renderer = vtkSmartPointer<vtkRenderer>::New();

    // VTK set background color
    this->m_renderer->SetBackground(0, 0, 0);

    // VTK/Qt wedded
    this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(m_renderer);

    // Set up action signals and slots
    connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
    connect(this->ui->actionSolve, SIGNAL(triggered()), this, SLOT(slotSolve()));

};

FemFlow::~FemFlow()
{
    // The smart pointers should clean up for up


}

void FemFlow::init()
{

    this->m_renderer->GetActiveCamera()->SetPosition(0.5, 0.5, 2);
    this->m_renderer->GetActiveCamera()->SetFocalPoint(0.5, 0.5, 0);

    /*
    vtkSmartPointer<vtkAxesActor> anAxesActor = vtkSmartPointer<vtkAxesActor>::New();

    vtkSmartPointer<vtkOrientationMarkerWidget> aWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    aWidget->SetOutlineColor(0.9300, 0.5700, 0.1300);
    aWidget->SetOrientationMarker(anAxesActor);
    aWidget->SetInteractor(this->ui->qvtkWidget->GetRenderWindow()->GetInteractor());
    aWidget->SetViewport(0.0, 0.0, 0.4, 0.4);
    aWidget->SetEnabled(1);
    aWidget->InteractiveOn();
    */

}

void FemFlow::slotExit()
{

    qApp->exit();

}

void FemFlow::slotSolve()
{

    this->solve1();
    //this->solve2();

}

void FemFlow::solve1()
{

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

    aBasicMesher.generate(aQuadrilateral, 40, 40, true);
    CMshMesh<double> aMesh = aBasicMesher.mesh();

    // Mesh visualization
    vtkSmartPointer<vtkPoints> sPoints = vtkSmartPointer<vtkPoints>::New();

    for (Integer i = 0; i < aMesh.nbNodes(); ++i)
    {

        unsigned int aNodeId = aMesh.nodeId(i);

        CMshNode<double> aNode = aMesh.node(aNodeId);

        sPoints->InsertPoint(aNodeId, aNode.x(), aNode.y(), aNode.z());

    }

    vtkSmartPointer<vtkUnstructuredGrid> anUnstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    anUnstructuredGrid->SetPoints(sPoints);

    for (Integer i = 0; i < aMesh.nbElements(); ++i)
    {

        unsigned int anElementId = aMesh.elementId(i);

        CMshElement<double> anElement = aMesh.element(anElementId);

        if (anElement.elementType() == ET_TRIANGLE)
        {
            vtkIdType ptIds[] = {
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(0))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(1))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(2)))
            };

            anUnstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, ptIds);
        }

    }

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

    double mu = 0.001;  // dynamic viscosity
    double rho = 1.0; // density

    double nu = mu / rho; // kinematic viscosity 

    std::cout << "Re = " << Ulid / nu << std::endl;

    SparseMatrix<double> G1 = gradient<double>(u, CP_X).matrixA;
    SparseMatrix<double> G2 = gradient<double>(v, CP_Y).matrixA;

    double dt = 5E-3;
    Integer nIter = 5000;

    this->m_renderer->RemoveActor(m_resultsActor);
    this->m_renderer->RemoveActor(m_scalarBarActor);

    this->m_resultsActor = vtkSmartPointer<vtkActor>::New();
    this->m_scalarBarActor = vtkSmartPointer<vtkScalarBarActor>::New();

    this->m_renderer->AddActor(m_resultsActor);
    this->m_renderer->AddActor(m_scalarBarActor);

    // Flow in a rectangle
    for (Integer i = 0; i < nIter; ++i)
    {

        std::cout << "Iteration = " << (i + 1) << std::endl;
        //std::cout << "Time = " << dt * (i + 1) << std::endl;

        // Velocity calculation
        CPdeEquation<double> aPdeEquationU(1.0 / dt * ddt<double>(u) -nu * laplacian<double>(u) + divergence<double>(u, v, dt) = 0);
        aPdeEquationU.setElimination(u);
        aPdeEquationU.solve(u);

        CPdeEquation<double> aPdeEquationV(1.0 / dt * ddt<double>(v) -nu * laplacian<double>(v) + divergence<double>(u, v, dt) = 0);
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

        vtkSmartPointer<vtkFloatArray> sScalars = vtkSmartPointer<vtkFloatArray>::New();

        sScalars->SetNumberOfTuples(aMesh.nbNodes());

        for (Integer i = 0; i < aMesh.nbNodes(); ++i)
        {

            unsigned int aNodeId = aMesh.nodeId(i);

            double vel = sqrt(u.u[aNodeId] * u.u[aNodeId] + v.u[aNodeId] * v.u[aNodeId]);

            sScalars->SetTuple1(aNodeId, vel);

        }

        anUnstructuredGrid->GetPointData()->SetScalars(sScalars);

        vtkSmartPointer<vtkGeometryFilter> aGeometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
        aGeometryFilter->SetInputData(anUnstructuredGrid);
        aGeometryFilter->Update();

        vtkSmartPointer<vtkBandedPolyDataContourFilter> aBoundedFilter = vtkSmartPointer<vtkBandedPolyDataContourFilter>::New();
        aBoundedFilter->SetInputData(aGeometryFilter->GetOutput());
        aBoundedFilter->GenerateValues(24, sScalars->GetRange());

        vtkSmartPointer<vtkLookupTable> aLookupTable = vtkSmartPointer<vtkLookupTable>::New();
        aLookupTable->SetNumberOfColors(256);
        aLookupTable->SetHueRange(0.667, 0.0);
        aLookupTable->Build();

        vtkSmartPointer<vtkPolyDataMapper> aBandedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        aBandedMapper->SetInputConnection(aBoundedFilter->GetOutputPort());
        aBandedMapper->SetScalarModeToUsePointData();
        aBandedMapper->SetScalarRange(sScalars->GetRange());
        aBandedMapper->SetLookupTable(aLookupTable);

        m_resultsActor->SetMapper(aBandedMapper);
        m_resultsActor->Modified();

        m_scalarBarActor->SetLookupTable(aLookupTable);
        m_scalarBarActor->SetTitle("Velocity");
        m_scalarBarActor->SetNumberOfLabels(4);
        m_scalarBarActor->Modified();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

        QCoreApplication::processEvents();

    }

}


void FemFlow::solve2()
{

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

    aBasicMesher.generate(aHexahedron, 25, 25, 1, true);
    CMshMesh<double> aMesh = aBasicMesher.mesh();

    // Mesh visualization
    vtkSmartPointer<vtkPoints> sPoints = vtkSmartPointer<vtkPoints>::New();

    for (Integer i = 0; i < aMesh.nbNodes(); ++i)
    {

        unsigned int aNodeId = aMesh.nodeId(i);

        CMshNode<double> aNode = aMesh.node(aNodeId);

        sPoints->InsertPoint(aNodeId, aNode.x(), aNode.y(), aNode.z());

    }

    vtkSmartPointer<vtkUnstructuredGrid> anUnstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    anUnstructuredGrid->SetPoints(sPoints);

    for (Integer i = 0; i < aMesh.nbElements(); ++i)
    {

        unsigned int anElementId = aMesh.elementId(i);

        CMshElement<double> anElement = aMesh.element(anElementId);

        if (anElement.elementType() == ET_TETRAHEDRON)
        {
            vtkIdType ptIds[] = {
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(0))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(1))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(2))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(3)))
            };

            anUnstructuredGrid->InsertNextCell(VTK_TETRA, 4, ptIds);
        }

    }

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

    double mu = 0.001;  // dynamic viscosity
    double rho = 1.0; // density

    double nu = mu / rho; // kinematic viscosity

    std::cout << "Re = " << Ulid / nu << std::endl;

    SparseMatrix<double> G1 = gradient<double>(u, CP_X).matrixA;
    SparseMatrix<double> G2 = gradient<double>(v, CP_Y).matrixA;

    double dt = 1E-2;
    Integer nIter = 5000;

    this->m_renderer->RemoveActor(m_resultsActor);
    this->m_renderer->RemoveActor(m_scalarBarActor);

    this->m_resultsActor = vtkSmartPointer<vtkActor>::New();
    this->m_scalarBarActor = vtkSmartPointer<vtkScalarBarActor>::New();

    this->m_renderer->AddActor(m_resultsActor);
    this->m_renderer->AddActor(m_scalarBarActor);

    // Flow in a box
    for (Integer i = 0; i < nIter; ++i)
    {

        std::cout << "Iteration = " << (i + 1) << std::endl;
        //std::cout << "Time = " << dt * (i + 1) << std::endl;

        // Velocity calculation
        CPdeEquation<double> aPdeEquationU(1.0 / dt * ddt<double>(u) -nu * laplacian<double>(u) +divergence<double>(u, v, dt) = 0);
        aPdeEquationU.setElimination(u);
        aPdeEquationU.solve(u);

        CPdeEquation<double> aPdeEquationV(1.0 / dt * ddt<double>(v) -nu * laplacian<double>(v) +divergence<double>(u, v, dt) = 0);
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

        vtkSmartPointer<vtkFloatArray> sScalars = vtkSmartPointer<vtkFloatArray>::New();

        sScalars->SetNumberOfTuples(aMesh.nbNodes());

        for (Integer i = 0; i < aMesh.nbNodes(); ++i)
        {

            unsigned int aNodeId = aMesh.nodeId(i);

            double vel = sqrt(u.u[aNodeId] * u.u[aNodeId] + v.u[aNodeId] * v.u[aNodeId]);

            sScalars->SetTuple1(aNodeId, vel);

        }

        anUnstructuredGrid->GetPointData()->SetScalars(sScalars);

        vtkSmartPointer<vtkPlane> aPlane = vtkSmartPointer<vtkPlane>::New();
        aPlane->SetOrigin(0, 0, 0);
        aPlane->SetNormal(0, 0, 1);

        vtkSmartPointer<vtkCutter> aCutter = vtkSmartPointer<vtkCutter>::New();
        aCutter->SetCutFunction(aPlane);
        aCutter->SetInputData(anUnstructuredGrid);
        aCutter->Update();

        vtkSmartPointer<vtkGeometryFilter> aGeometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
        aGeometryFilter->SetInputConnection(aCutter->GetOutputPort());
        aGeometryFilter->Update();

        vtkSmartPointer<vtkBandedPolyDataContourFilter> aBoundedFilter = vtkSmartPointer<vtkBandedPolyDataContourFilter>::New();
        aBoundedFilter->SetInputData(aGeometryFilter->GetOutput());
        aBoundedFilter->GenerateValues(24, sScalars->GetRange());

        vtkSmartPointer<vtkLookupTable> aLookupTable = vtkSmartPointer<vtkLookupTable>::New();
        aLookupTable->SetNumberOfColors(256);
        aLookupTable->SetHueRange(0.667, 0.0);
        aLookupTable->Build();

        vtkSmartPointer<vtkPolyDataMapper> aBandedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        aBandedMapper->SetInputConnection(aBoundedFilter->GetOutputPort());
        aBandedMapper->SetScalarModeToUsePointData();
        aBandedMapper->SetScalarRange(sScalars->GetRange());
        aBandedMapper->SetLookupTable(aLookupTable);

        m_resultsActor->SetMapper(aBandedMapper);
        m_resultsActor->Modified();

        m_scalarBarActor->SetLookupTable(aLookupTable);
        m_scalarBarActor->SetTitle("Velocity");
        m_scalarBarActor->SetNumberOfLabels(4);
        m_scalarBarActor->Modified();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

        QCoreApplication::processEvents();

    }
}

