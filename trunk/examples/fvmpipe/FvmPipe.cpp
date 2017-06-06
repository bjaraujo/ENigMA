// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#include "ui_FvmPipe.h"
#include "FvmPipe.h"

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
FvmPipe::FvmPipe()
{

    this->ui = new Ui_FvmPipe;
    this->ui->setupUi(this);

    this->setWindowTitle("Fvm Pipe");

    // VTK Renderer
    this->m_renderer = vtkSmartPointer<vtkRenderer>::New();

    // Add Actor to renderer
    vtkSmartPointer<vtkAxesActor> anAxesActor = vtkSmartPointer<vtkAxesActor>::New();
    this->m_renderer->AddActor(anAxesActor);

    // VTK set background color
    this->m_renderer->SetBackground(0, 0, 0);

    // VTK/Qt wedded
    this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(m_renderer);

    // Set up action signals and slots
    connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
    connect(this->ui->actionSolve, SIGNAL(triggered()), this, SLOT(slotSolve()));

};

FvmPipe::~FvmPipe()
{
    // The smart pointers should clean up for up


}

void FvmPipe::init()
{

    this->m_renderer->GetActiveCamera()->SetPosition(1.0, 0.0, 0.25);
    this->m_renderer->GetActiveCamera()->SetFocalPoint(0.0, 0.0, 0.25);

}

void FvmPipe::slotExit()
{

    qApp->exit();

}

void FvmPipe::slotSolve()
{

    this->solve();

}

void FvmPipe::solve()
{

    CPosGmsh<double> aGmsPos;
    CPdeField<double> aField;

    aGmsPos.load(aField, "pipe.msh");

    CMshMesh<double> aMesh = aField.mesh();

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

    aMesh.generateFaces(1E-8);

    aMesh.calculateFaceCentroid();
    aMesh.calculateElementCentroid();

    CFvmMesh<double> aFvmMesh(aMesh);

    double mu = 1.0;  // dynamic viscosity
    double rho = 1.0; // density

    double nu = mu / rho; // kinematic viscosity 

    double U = 0.1;
    double d = 0.1;
    double L = 0.5;

    // Hagen Poiseuille Flow
    double dp = 32 * mu * L / (d * d) * U;

    std::cout << "dp = " << dp << std::endl;
    std::cout << "Umax = " << (d * d) * dp / (32 * mu * L) * 2 << std::endl;
    
    CFvmPisoSolver<double> aPisoSolver(aFvmMesh);

    aPisoSolver.setGravity(0.0, 0.0, 0.0);

    aPisoSolver.setMaterialProperties(rho, mu);

    std::vector<Integer> sWallFaceIds;

    sWallFaceIds.clear();
    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<double> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (aFace.centroid().z() > 0.0 + 1E-6 && aFace.centroid().z() < 0.5 - 1E-6)
        {
            sWallFaceIds.push_back(aFaceId);
        }

    }

    aPisoSolver.setBoundaryVelocity(sWallFaceIds, BT_WALL_NO_SLIP, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(sWallFaceIds, BT_WALL_NO_SLIP, 0.0);

    std::vector<Integer> sInletFaceIds;

    sInletFaceIds.clear();
    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<double> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (aFace.centroid().z() < 0.0 + 1E-6)
        {
            sInletFaceIds.push_back(aFaceId);
        }

    }

    aPisoSolver.setBoundaryVelocity(sInletFaceIds, BT_INLET_FLOW, 0.0, 0.0, U);
    aPisoSolver.setBoundaryPressure(sInletFaceIds, BT_INLET_FLOW, 0.0);

    //aPisoSolver.setBoundaryPressure(sInletFaceIds, BT_INLET_PRESSURE, dp);

    std::vector<Integer> sOutletFaceIds;

    sOutletFaceIds.clear();
    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<double> aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (aFace.centroid().z() > 0.5 - 1E-6)
        {
            sOutletFaceIds.push_back(aFaceId);
        }

    }

    aPisoSolver.setBoundaryVelocity(sOutletFaceIds, BT_OUTLET, 0.0, 0.0, 0.0);
    aPisoSolver.setBoundaryPressure(sOutletFaceIds, BT_OUTLET, 0.0);

    this->m_renderer->RemoveActor(m_resultsActor);
    this->m_renderer->RemoveActor(m_scalarBarActor);

    this->m_resultsActor = vtkSmartPointer<vtkActor>::New();
    this->m_scalarBarActor = vtkSmartPointer<vtkScalarBarActor>::New();

    this->m_renderer->AddActor(m_resultsActor);
    this->m_renderer->AddActor(m_scalarBarActor);

    double dt = 1E-5;
    Integer nIter = 500;
    Integer nInitSteps = 10;

    // Flow in a pipe
    for (Integer i = 0; i < nIter; ++i)
    {

        if (i < nInitSteps)
        {
            std::cout << "Iteration = " << (i + 1) << std::endl;
            aPisoSolver.iterate(dt, true);
        }
        else
        {
            std::cout << "Time = " << dt * (i - nInitSteps + 1) << std::endl;
            aPisoSolver.iterate(dt, false);
        }

        double massError;
        aPisoSolver.checkMassConservation(massError);
        std::cout << "Mass conservation error: " << std::scientific << massError << std::endl;

        vtkSmartPointer<vtkFloatArray> sScalars = vtkSmartPointer<vtkFloatArray>::New();

        sScalars->SetNumberOfTuples(aFvmMesh.nbControlVolumes());

        double p_max = std::numeric_limits<double>::min();

        for (Integer i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
        {

            unsigned int aControlVolumeId = aFvmMesh.controlVolumeId(i);

            double u = aPisoSolver.u(aControlVolumeId);
            double v = aPisoSolver.v(aControlVolumeId);
            double w = aPisoSolver.w(aControlVolumeId);

            double vel = sqrt(u*u + v*v + w*w);

            p_max = std::max(p_max, aPisoSolver.p(aControlVolumeId));

            //CGeoVector<double> gradp = aPisoSolver.gradp(aControlVolumeId);
            //double grad = sqrt(gradp.x()*gradp.x() + gradp.y()*gradp.y() + gradp.z()*gradp.z());

            //sScalars->SetTuple1(aControlVolumeId, vel);
            //sScalars->SetTuple1(aControlVolumeId, aPisoSolver.u(aControlVolumeId));
            //sScalars->SetTuple1(aControlVolumeId, aPisoSolver.v(aControlVolumeId));
            //sScalars->SetTuple1(aControlVolumeId, aPisoSolver.w(aControlVolumeId));
            sScalars->SetTuple1(aControlVolumeId, aPisoSolver.p(aControlVolumeId));
            //sScalars->SetTuple1(aControlVolumeId, grad);
            //sScalars->SetTuple1(aControlVolumeId, aPisoSolver.massError(aControlVolumeId));

        }

        std::cout << "pmax = " << std::fixed << p_max << std::endl;

        double Ui_max = std::numeric_limits<double>::min();        
        for (Integer i = 0; i < sInletFaceIds.size(); i++)
        {
            Ui_max = std::max(Ui_max, aPisoSolver.wf(sInletFaceIds[i]));
        }
        std::cout << "Inlet Umax = " << std::fixed << Ui_max << std::endl;

        double Uo_max = std::numeric_limits<double>::min();
        for (Integer i = 0; i < sOutletFaceIds.size(); i++)
        {
            Uo_max = std::max(Uo_max, aPisoSolver.wf(sOutletFaceIds[i]));
        }
        std::cout << "Outlet Umax = " << std::fixed << Uo_max << std::endl;

        anUnstructuredGrid->GetCellData()->SetScalars(sScalars);

        vtkSmartPointer<vtkPlane> aPlane = vtkSmartPointer<vtkPlane>::New();
        aPlane->SetOrigin(0, 0, 0);
        aPlane->SetNormal(1, 0, 0);

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
        //m_scalarBarActor->SetTitle("Mass Error");
        //m_scalarBarActor->SetTitle("Velocity");
        m_scalarBarActor->SetTitle("Pressure");
        m_scalarBarActor->SetNumberOfLabels(4);
        m_scalarBarActor->Modified();

        this->ui->qvtkWidget->GetRenderWindow()->Render();
        
        QCoreApplication::processEvents();
        
    }

}

