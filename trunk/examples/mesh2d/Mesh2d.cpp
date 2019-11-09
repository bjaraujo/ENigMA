// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "ui_Mesh2d.h"
#include "Mesh2d.h"

#include <QElapsedTimer>

#include <vtkEventQtSlotConnect.h>
#include <vtkDataObjectToTable.h>
#include <vtkQtTableView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// Constructor
Mesh2D::Mesh2D()
{

    this->ui = new Ui_Mesh2d;
    this->ui->setupUi(this);

    // Qt Table View
    this->m_tableView = vtkSmartPointer<vtkQtTableView>::New();

    // Place the table view in the designer form
    this->ui->tableFrame->layout()->addWidget(this->m_tableView->GetWidget());

    // Mapper
    VTK_CREATE(vtkDataSetMapper, mapper);

    VTK_CREATE(vtkUnstructuredGrid, unstructuredGrid);
    mapper->SetInputData(unstructuredGrid);

    mapper->ImmediateModeRenderingOn();

    // Actor in scene
    VTK_CREATE(vtkActor, actor);
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    actor->GetProperty()->LightingOn();

    VTK_CREATE(vtkAxesActor, axes);

    // VTK Renderer
    VTK_CREATE(vtkRenderer, ren);

    // Add Actor to renderer
    ren->AddActor(axes);
    ren->AddActor(actor);

    // VTK/Qt wedded
    this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(ren);

    // Just a bit of Qt interest: Culling off the
    // point data and handing it to a vtkQtTableView
    //VTK_CREATE(vtkDataObjectToTable, toTable);
    //toTable->SetInputConnection(elevation->GetOutputPort());
    //toTable->SetFieldType(vtkDataObjectToTable::POINT_DATA);

    // Here we take the end of the VTK pipeline and give it to a Qt View
    //this->m_tableView->SetRepresentationFromInputConnection(toTable->GetOutputPort());

    // Set up action signals and slots
    connect(this->ui->actionOpenFile, SIGNAL(triggered()), this, SLOT(slotOpenFile()));
    connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));

    m_connections = vtkEventQtSlotConnect::New();
    this->m_connections->Connect(this->ui->qvtkWidget->GetRenderWindow()->GetInteractor(),
        vtkCommand::KeyPressEvent, this,
        SLOT(slotKeyPressed(vtkObject*, unsigned long, void*, void*)), 0, 1.0);

    m_mesher.setIntervals(9999999, 9999999);

    m_mesher.onUpdate = std::bind1st(std::mem_fun(&Mesh2D::updateData), this);

};

Mesh2D::~Mesh2D()
{
    // The smart pointers should clean up for up


}

void Mesh2D::slotKeyPressed(vtkObject *, unsigned long, void *, void *)
{

    vtkRenderWindowInteractor *rwi = this->ui->qvtkWidget->GetRenderWindow()->GetInteractor();

    std::string key = rwi->GetKeySym();

    // Output the key that was pressed

    //std::cout << "Pressed " << key << std::endl;

    if (key == "Right")
    {

        drawMesh();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    }
    else if (key == "Left")
    {

        drawMesh();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    }
    else if (key == "m" || key == "M")
    {

        generateMesh(60, 30);
        drawMesh();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    }
    else if (key == "w" || key == "W")
    {

        vtkRenderer* ren = this->ui->qvtkWidget->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

        vtkActor* actor = ren->GetActors()->GetLastActor();

        actor->GetProperty()->SetRepresentationToWireframe();
        actor->GetProperty()->LightingOff();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    }
    else if (key == "s" || key == "S")
    {

        vtkRenderer* ren = this->ui->qvtkWidget->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

        vtkActor* actor = ren->GetActors()->GetLastActor();

        actor->GetProperty()->SetRepresentationToSurface();
        actor->GetProperty()->LightingOn();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    }

}

// Action to be taken upon file open 
void Mesh2D::slotOpenFile()
{

}

void Mesh2D::slotExit()
{

    qApp->exit();

}

bool Mesh2D::generateMesh(const unsigned int nu, const unsigned int nv)
{

    std::cout << "Start..." << std::endl;

    const double d = 0.25;

    std::cout << "nu = " << nu << std::endl;
    std::cout << "nv = " << nv << std::endl;

    std::cout << "d = " << d << std::endl;

    CGeoCoordinate<double> aVertex1(0.0*nu*d, 0.0*nv*d, 0.0);
    CGeoCoordinate<double> aVertex2(1.0*nu*d, 0.0*nv*d, 0.0);
    CGeoCoordinate<double> aVertex3(1.0*nu*d, 1.0*nv*d, 0.0);
    CGeoCoordinate<double> aVertex4(0.0*nu*d, 1.0*nv*d, 0.0);

    CGeoCoordinate<double> aVertex5(0.1*nu*d, 0.1*nv*d, 0.0);
    CGeoCoordinate<double> aVertex6(0.5*nu*d, 0.1*nv*d, 0.0);
    CGeoCoordinate<double> aVertex7(0.5*nu*d, 0.5*nv*d, 0.0);
    CGeoCoordinate<double> aVertex8(0.1*nu*d, 0.5*nv*d, 0.0);

    CGeoQuadrilateral<double> aQuadrilateral1;
    CGeoQuadrilateral<double> aQuadrilateral2;

    aQuadrilateral1.addVertex(aVertex1);
    aQuadrilateral1.addVertex(aVertex2);
    aQuadrilateral1.addVertex(aVertex3);
    aQuadrilateral1.addVertex(aVertex4);

    aQuadrilateral2.addVertex(aVertex8);
    aQuadrilateral2.addVertex(aVertex7);
    aQuadrilateral2.addVertex(aVertex6);
    aQuadrilateral2.addVertex(aVertex5);

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aQuadrilateral1, nu / 4, nv / 4, false);
    CMshMesh<double> aMesh1 = aBasicMesher.mesh();

    aBasicMesher.generate(aQuadrilateral2, nu / 4, nv / 4, false);
    CMshMesh<double> aMesh2 = aBasicMesher.mesh();

    CMshMesh<double> aBoundaryMesh1 = aMesh1.extractBoundary(1E-6);
    CMshMesh<double> aBoundaryMesh2 = aMesh2.extractBoundary(1E-6);

    CMshMesh<double> aBoundaryMesh3;

    aBoundaryMesh3.addMesh(aBoundaryMesh1);
    aBoundaryMesh3.addMesh(aBoundaryMesh2);

    aBoundaryMesh3.generateFaces(1E-4);

    std::vector<CGeoCoordinate<double> > sInteriorPoints;

    bool res = false;

    try
    {

        QElapsedTimer timer;

        timer.start();

        res = m_mesher.generate(aBoundaryMesh3, 9999, sInteriorPoints, d, 0.05, 1E-1);

        std::cout << "Done. Elapsed time = " << timer.elapsed() / 1000.0 << " s" << std::endl;
        std::cout << "Number of elements: " << m_mesher.mesh().nbElements() << std::endl;

        /*
        for (int i = 0; i < 5; i++)
        {
            m_mesher.flipEdges(1E-5);
            m_mesher.relaxNodes(1E-5);
        }
        */

    }
    catch (const std::exception& e)
    {
        std::cout << "Error: std exception: " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cout << "Error: unknown exception" << std::endl;
    }

    CPosGmsh<double> aPosGmsh;
    CPdeField<double> aField;

    aField.setMesh(m_mesher.mesh());
    aPosGmsh.save(aField, "tri_surface.msh", "tris");

    return res;

}

void Mesh2D::drawMesh()
{

    if (m_mesher.mesh().nbNodes() == 0)
        return;

    VTK_CREATE(vtkPoints, points);

    for (Integer i = 0; i < m_mesher.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = m_mesher.mesh().nodeId(i);

        CMshNode<double> aNode = m_mesher.mesh().node(aNodeId);

        points->InsertNextPoint(aNode.x(), aNode.y(), aNode.z());

    }

    VTK_CREATE(vtkUnstructuredGrid, unstructuredGrid);
    unstructuredGrid->SetPoints(points);

    for (Integer i = 0; i < m_mesher.mesh().nbElements(); ++i)
    {

        Integer anElementId = m_mesher.mesh().elementId(i);

        CMshElement<double> anElement = m_mesher.mesh().element(anElementId);

        if (anElement.elementType() == ET_TRIANGLE)
        {
            vtkIdType ptIds[] = { static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(0))),
                static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(1))),
                static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(2))) };

            unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, ptIds);
        }

        if (anElement.elementType() == ET_QUADRILATERAL)
        {
            vtkIdType ptIds[] = { static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(0))),
                static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(1))),
                static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(2))),
                static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(3))) };

            unstructuredGrid->InsertNextCell(VTK_QUAD, 4, ptIds);
        }

    }

    vtkRenderer* ren = this->ui->qvtkWidget->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

    vtkMapper* mapper = ren->GetActors()->GetLastActor()->GetMapper();

    (reinterpret_cast<vtkDataSetMapper*>(mapper))->SetInputData(unstructuredGrid);

}

int Mesh2D::updateData(int)
{

    m_mutex.lock();

    drawMesh();

    this->ui->qvtkWidget->GetRenderWindow()->Render();

    QCoreApplication::processEvents();

    m_mutex.unlock();

    return 0;

}

