// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "ui_FvmFlow.h"
#include "FvmFlow.h"

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

#include "MshBasicMesher.hpp"
#include "SleSystem.hpp"
#include "PdeEquation.hpp"
#include "PosGnuplot.hpp"
#include "PosGmsh.hpp"
#include "AnaTemperature.hpp"
#include "FvmPisoSolver.hpp"

using namespace ENigMA::fvm;
using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::analytical;

// Constructor
FvmFlow::FvmFlow()
{

    this->ui = new Ui_FvmFlow;
    this->ui->setupUi(this);

    this->setWindowTitle("Fvm Flow");

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

FvmFlow::~FvmFlow()
{
    // The smart pointers should clean up for up


}

void FvmFlow::init()
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

void FvmFlow::slotExit()
{

    qApp->exit();

}

void FvmFlow::slotSolve()
{

    this->solve();

}

void FvmFlow::solve()
{

    CGeoCoordinate<double> aVertex1(+0.00, +0.00, -5.0);
    CGeoCoordinate<double> aVertex2(+1.00, +0.00, -5.0);
    CGeoCoordinate<double> aVertex3(+1.00, +1.00, -5.0);
    CGeoCoordinate<double> aVertex4(+0.00, +1.00, -5.0);
    CGeoCoordinate<double> aVertex5(+0.00, +0.00, +5.0);
    CGeoCoordinate<double> aVertex6(+1.00, +0.00, +5.0);
    CGeoCoordinate<double> aVertex7(+1.00, +1.00, +5.0);
    CGeoCoordinate<double> aVertex8(+0.00, +1.00, +5.0);

    CGeoHexahedron<double> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);
    
    Integer nx = 80;
    Integer ny = 80;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nx, ny, 1);
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

        if (anElement.elementType() == ET_HEXAHEDRON)
        {
            vtkIdType ptIds[] = {
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(3))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(2))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(1))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(0))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(7))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(6))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(5))),
                static_cast<vtkIdType> (aMesh.nodeIndex(anElement.nodeId(4)))
            };

            anUnstructuredGrid->InsertNextCell(VTK_HEXAHEDRON, 8, ptIds);
        }

    }

    aMesh.generateFaces(1E-12);

    aMesh.calculateFaceCentroid();
    aMesh.calculateElementCentroid();

    CFvmMesh<double> aFvmMesh(aMesh);

    double U = 1.0; // Lid velocity
    double mu = 0.001; // dynamic viscosity
    double rho = 1.0; // density

    double nu = mu / rho; // kinematic viscosity 

    std::cout << "Re = " << U / nu << std::endl;

    CFvmPisoSolver<double> aPisoSolver(aFvmMesh);

    aPisoSolver.setGravity(0.0, 0.0, 0.0);
    aPisoSolver.setMaterialProperties(rho, mu);

    std::vector<Integer> sFaceIds;

    sFaceIds.clear();
    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<double> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (aFace.centroid().x() == 0.0 ||
            aFace.centroid().x() == 1.0 ||
            aFace.centroid().y() == 0.0 ||
            aFace.centroid().y() == 1.0)
        {
            sFaceIds.push_back(aFaceId);
        }

    }

    aPisoSolver.setBoundaryVelocity(sFaceIds, BT_WALL_NO_SLIP, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(sFaceIds, BT_WALL_NO_SLIP, 0.0);

    sFaceIds.clear();
    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<double> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (aFace.centroid().y() == 1.0)
        {
            sFaceIds.push_back(aFaceId);
        }

    }

    aPisoSolver.setBoundaryVelocity(sFaceIds, BT_WALL_NO_SLIP, U, 0.0, 0.0);

    double dt = U / nx;            // Courant < 1
    Integer nIter = 5000;

    this->m_renderer->RemoveActor(m_resultsActor);
    this->m_renderer->RemoveActor(m_scalarBarActor);

    this->m_resultsActor = vtkSmartPointer<vtkActor>::New();
    this->m_scalarBarActor = vtkSmartPointer<vtkScalarBarActor>::New();

    this->m_renderer->AddActor(m_resultsActor);
    this->m_renderer->AddActor(m_scalarBarActor);

    double r0, r = 0.0;

    // Flow in a rectangle
    for (Integer ii = 0; ii < nIter; ++ii)
    {

        std::cout << "Iteration = " << (ii + 1) << std::endl;
        std::cout << "Time = " << dt * (ii + 1) << std::endl;

        aPisoSolver.iterate(dt);

        double massError;
        aPisoSolver.checkMassConservation(massError);
        std::cout << "Mass conservation error: " << std::scientific << massError << std::endl;

        vtkSmartPointer<vtkFloatArray> sScalars = vtkSmartPointer<vtkFloatArray>::New();

        sScalars->SetNumberOfTuples(aFvmMesh.nbControlVolumes());

        for (Integer i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
        {

            unsigned int aControlVolumeId = aFvmMesh.controlVolumeId(i);

            double u = aPisoSolver.u(aControlVolumeId);
            double v = aPisoSolver.v(aControlVolumeId);

            double vel = sqrt(u*u + v*v);

            sScalars->SetTuple1(aControlVolumeId, vel);

        }

        anUnstructuredGrid->GetCellData()->SetScalars(sScalars);

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

        vtkSmartPointer<vtkCellDataToPointData> aCellDataToPointDataFilter1 = vtkSmartPointer<vtkCellDataToPointData>::New();
        aCellDataToPointDataFilter1->SetInputConnection(aGeometryFilter->GetOutputPort());
        aCellDataToPointDataFilter1->Update();

        vtkSmartPointer<vtkBandedPolyDataContourFilter> aBoundedFilter = vtkSmartPointer<vtkBandedPolyDataContourFilter>::New();
        aBoundedFilter->SetInputData(aCellDataToPointDataFilter1->GetOutput());
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

        double ru, rv, rw, rp;
        aPisoSolver.residual(ru, rv, rw, rp);

        r0 = r;
        r = sqrt(ru * ru + rv * rv);

        std::cout << "Residual: " << std::scientific << r << std::endl;

    }

}

