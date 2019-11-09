// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "ui_FvmCicsam.h"
#include "FvmCicsam.h"

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
#include "FvmVofSolver.hpp"

using namespace ENigMA::fvm;
using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::analytical;

// Constructor
FvmCicsam::FvmCicsam()
{

    this->ui = new Ui_FvmCicsam;
    this->ui->setupUi(this);

    this->setWindowTitle("Fvm CICSAM");

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

FvmCicsam::~FvmCicsam()
{
    // The smart pointers should clean up for up


}

void FvmCicsam::init()
{

    this->m_renderer->GetActiveCamera()->SetPosition(0.5, 0.5, 2);
    this->m_renderer->GetActiveCamera()->SetFocalPoint(0.5, 0.5, 0);

}

void FvmCicsam::slotExit()
{

    qApp->exit();

}

void FvmCicsam::slotSolve()
{

    this->solve();

}

void FvmCicsam::solve()
{

    CGeoCoordinate<double> aVertex1(+0.00, +0.00, -0.05);
    CGeoCoordinate<double> aVertex2(+1.00, +0.00, -0.05);
    CGeoCoordinate<double> aVertex3(+1.00, +0.10, -0.05);
    CGeoCoordinate<double> aVertex4(+0.00, +0.10, -0.05);
    CGeoCoordinate<double> aVertex5(+0.00, +0.00, +0.05);
    CGeoCoordinate<double> aVertex6(+1.00, +0.00, +0.05);
    CGeoCoordinate<double> aVertex7(+1.00, +0.10, +0.05);
    CGeoCoordinate<double> aVertex8(+0.00, +0.10, +0.05);

    CGeoHexahedron<double> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    Integer nx = 50;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aHexahedron, nx, 1, 1);
    CMshMesh<double> aOriginalMesh = aBasicMesher.mesh();

    aOriginalMesh.generateFaces(1E-4);

    aOriginalMesh.calculateFaceCentroid();
    aOriginalMesh.calculateElementCentroid();

    // Apply gamma
    CFvmMesh<double> aFvmMesh(aOriginalMesh);

    CFvmVofSolver<double> aVofSolver(aFvmMesh);

    for (Integer i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
    {

        Integer aControlVolumeId = aFvmMesh.controlVolumeId(i);

        CFvmControlVolume<double>& aControlVolume = aFvmMesh.controlVolume(aControlVolumeId);

        aControlVolume.calculateCentroid();

        if (aControlVolume.centroid().x() > 0.4 && aControlVolume.centroid().x() < 0.6)
            aVofSolver.setInitialGamma(aControlVolumeId, 1.0);
        else
            aVofSolver.setInitialGamma(aControlVolumeId, 0.0);
        
    }

    // Mesh visualization
    vtkSmartPointer<vtkPoints> sPoints = vtkSmartPointer<vtkPoints>::New();

    for (Integer i = 0; i < aOriginalMesh.nbNodes(); ++i)
    {

        Integer aNodeId = aOriginalMesh.nodeId(i);

        CMshNode<double>& aNode = aOriginalMesh.node(aNodeId);

        sPoints->InsertPoint(aNodeId, aNode.x(), aNode.y(), aNode.z());

    }

    vtkSmartPointer<vtkUnstructuredGrid> anUnstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    anUnstructuredGrid->SetPoints(sPoints);

    for (Integer i = 0; i < aOriginalMesh.nbElements(); ++i)
    {

        Integer anElementId = aOriginalMesh.elementId(i);

        CMshElement<double>& anElement = aOriginalMesh.element(anElementId);

        if (anElement.elementType() == ET_HEXAHEDRON)
        {
            vtkIdType ptIds[] = {
                static_cast<vtkIdType> (aOriginalMesh.nodeIndex(anElement.nodeId(3))),
                static_cast<vtkIdType> (aOriginalMesh.nodeIndex(anElement.nodeId(2))),
                static_cast<vtkIdType> (aOriginalMesh.nodeIndex(anElement.nodeId(1))),
                static_cast<vtkIdType> (aOriginalMesh.nodeIndex(anElement.nodeId(0))),
                static_cast<vtkIdType> (aOriginalMesh.nodeIndex(anElement.nodeId(7))),
                static_cast<vtkIdType> (aOriginalMesh.nodeIndex(anElement.nodeId(6))),
                static_cast<vtkIdType> (aOriginalMesh.nodeIndex(anElement.nodeId(5))),
                static_cast<vtkIdType> (aOriginalMesh.nodeIndex(anElement.nodeId(4)))
            };

            anUnstructuredGrid->InsertNextCell(VTK_HEXAHEDRON, 8, ptIds);
        }

    }

    double rho0 = 1.0; // density
    double mu0 = 1.0; // dynamic viscosity

    double rho1 = 1.0; // density
    double mu1 = 1.0; // dynamic viscosity

    aVofSolver.setGravity(0.0, 0.0, 0.0);

    aVofSolver.setMaterialProperties(rho0, mu0, rho1, mu1);

    std::vector<Integer> sInletIds;

    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<double>& aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (aFace.centroid().x() < 1E-6)
        {
            sInletIds.push_back(aFaceId);
        }

    }

    aVofSolver.setBoundaryVelocity(sInletIds, BT_INLET_FLOW, 1.0, 0.0, 0.0);
    aVofSolver.setBoundaryPressure(sInletIds, BT_INLET_FLOW, 0.0);
    aVofSolver.setBoundaryGamma(sInletIds, BT_INLET_FLOW, 0.0);

    std::vector<Integer> sOutletIds;

    for (Integer i = 0; i < aFvmMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aFvmMesh.faceId(i);

        CFvmFace<double>& aFace = aFvmMesh.face(aFaceId);

        aFace.calculateCentroid();

        if (aFace.centroid().x() > (1.0 - 1E-6))
        {
            sOutletIds.push_back(aFaceId);
        }

    }

    aVofSolver.setBoundaryVelocity(sOutletIds, BT_OUTLET, 0.0, 0.0, 0.0);
    aVofSolver.setBoundaryPressure(sOutletIds, BT_OUTLET, 0.0);
    aVofSolver.setBoundaryGamma(sOutletIds, BT_OUTLET, 0.0);

    double dt = 1E-4;
    Integer nIter = 1000;

    this->m_renderer->RemoveActor(m_resultsActor);
    this->m_renderer->RemoveActor(m_scalarBarActor);

    this->m_resultsActor = vtkSmartPointer<vtkActor>::New();
    this->m_scalarBarActor = vtkSmartPointer<vtkScalarBarActor>::New();

    this->m_renderer->AddActor(m_resultsActor);
    this->m_renderer->AddActor(m_scalarBarActor);

    // Flow in a rectangle
    for (Integer ii = 0; ii < nIter; ++ii)
    {

        //std::cout << "Iteration = " << (ii + 1) << std::endl;
        std::cout << "Time = " << std::scientific << dt * (ii + 1) << std::endl;
        //std::cout << "Interval = " << dt << std::endl;

        aVofSolver.iterate(dt);

        double massError;
        aVofSolver.checkMassConservation(massError);
        std::cout << "Mass conservation error: " << std::scientific << massError << std::endl;

        vtkSmartPointer<vtkFloatArray> sScalars = vtkSmartPointer<vtkFloatArray>::New();

        sScalars->SetNumberOfTuples(aFvmMesh.nbControlVolumes());

        for (Integer i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
        {

            Integer aControlVolumeId = aFvmMesh.controlVolumeId(i);

            sScalars->SetTuple1(aControlVolumeId, aVofSolver.s(aControlVolumeId));

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
        m_scalarBarActor->SetTitle("Gamma");
        m_scalarBarActor->SetNumberOfLabels(4);
        m_scalarBarActor->Modified();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

        QCoreApplication::processEvents();

    }

}

