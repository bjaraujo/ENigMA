// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#include "ui_StlUtils.h"
#include "StlUtils.h"

#include <vtkEventQtSlotConnect.h>
#include <vtkDataObjectToTable.h>
#include <vtkQtTableView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkSTLReader.h>
#include <vtkFeatureEdges.h>

// Constructor
StlUtils::StlUtils() 
{
    this->ui = new Ui_StlUtils;
    this->ui->setupUi(this);

    this->setWindowTitle("Stl Reader");

    this->m_stlReader = vtkSmartPointer<vtkSTLReader>::New();

    // Qt Table View
    this->m_tableView = vtkSmartPointer<vtkQtTableView>::New();

    // Place the table view in the designer form
    this->ui->tableFrame->layout()->addWidget(this->m_tableView->GetWidget());

    // Mappers
    this->m_stlMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    this->m_stlMapper->ImmediateModeRenderingOn();

    // Actor in scene
    this->m_axesActor = vtkSmartPointer<vtkAxesActor>::New();

    this->m_stlActor = vtkSmartPointer<vtkActor>::New();
    this->m_stlActor->SetMapper(m_stlMapper);
    this->m_stlActor->GetProperty()->SetRepresentationToSurface();
    this->m_stlActor->GetProperty()->LightingOn();

    this->m_edgeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    this->m_edgeMapper->ImmediateModeRenderingOn();

    this->m_edgeActor = vtkSmartPointer<vtkActor>::New();
    this->m_edgeActor->SetMapper(m_edgeMapper);

    this->m_pointMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    this->m_pointMapper->ImmediateModeRenderingOn();

    this->m_pointActor = vtkSmartPointer<vtkActor>::New();
    this->m_pointActor->SetMapper(m_pointMapper);
    this->m_pointActor->GetProperty()->SetPointSize(5);

    // VTK Renderer
    this->m_renderer = vtkSmartPointer<vtkRenderer>::New();

    // Add Actor to renderer
    this->m_renderer->AddActor(this->m_axesActor);
    this->m_renderer->AddActor(this->m_stlActor);
    this->m_renderer->AddActor(this->m_edgeActor);
    this->m_renderer->AddActor(this->m_pointActor);

    // VTK/Qt wedded
    this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(m_renderer);

    // Set up action signals and slots
    connect(this->ui->actionOpenFile, SIGNAL(triggered()), this, SLOT(slotOpenFile()));
    connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));

    this->m_connections = vtkEventQtSlotConnect::New();
    this->m_connections->Connect(this->ui->qvtkWidget->GetRenderWindow()->GetInteractor(), 
        vtkCommand::KeyPressEvent, this, 
        SLOT(slotKeyPressed(vtkObject*, unsigned long, void*, void*)), 0, 1.0);

};

StlUtils::~StlUtils()
{
    // The smart pointers should clean up for up


}

void StlUtils::slotKeyPressed(vtkObject *, unsigned long, void *, void *)
{

    vtkRenderWindowInteractor *rwi = this->ui->qvtkWidget->GetRenderWindow()->GetInteractor();

    std::string key = rwi->GetKeySym();

}

// Action to be taken upon file open 
void StlUtils::slotOpenFile()
{

    QString selfilter = tr("STL files (*.stl)");

    m_stlFileName = QFileDialog::getOpenFileName(this, tr("Import"), "", tr("STL files (*.stl)"), &selfilter);

    if (m_stlFileName == QString::null) return;

    QFileInfo fileInfo(m_stlFileName);

    if (!QFileInfo(m_stlFileName).exists())
    {
        qDebug() << "Error: file does not exist: " << m_stlFileName;
        return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);

    QString fileExt = QFileInfo(QString(m_stlFileName)).completeSuffix();

    if (QString::compare(fileExt, "stl", Qt::CaseInsensitive) == 0)
    {

        this->readStlFile();

        this->drawStlFile();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    }

    QApplication::restoreOverrideCursor();

}

void StlUtils::slotExit() 
{

    qApp->exit();

}

void StlUtils::readStlFile()
{

    this->m_stlReader->SetFileName(m_stlFileName.toStdString().c_str());
    this->m_stlReader->Update();

    qDebug() << "Number of facets: " << m_stlReader->GetOutput()->GetNumberOfPolys();

}

void StlUtils::drawStlFile()
{

    this->m_stlMapper->SetInputConnection(m_stlReader->GetOutputPort());

    this->m_boundaryEdges = vtkSmartPointer<vtkFeatureEdges>::New();
    this->m_boundaryEdges->SetInputConnection(m_stlReader->GetOutputPort());
    this->m_boundaryEdges->BoundaryEdgesOn();
    this->m_boundaryEdges->FeatureEdgesOff();
    this->m_boundaryEdges->ManifoldEdgesOff();
    this->m_boundaryEdges->NonManifoldEdgesOff();
    this->m_boundaryEdges->Update();

    this->m_edgeMapper->SetInputConnection(m_boundaryEdges->GetOutputPort());

    this->m_renderer->ResetCamera();

}
