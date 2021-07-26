
#include "ui_Mesh3d.h"
#include "Mesh3d.h"

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
Mesh3D::Mesh3D(QString fileName) 
{

    m_fileName = fileName;

    this->ui = new Ui_Mesh3d;
    this->ui->setupUi(this);

    // Qt Table View
    this->m_tableView = vtkSmartPointer<vtkQtTableView>::New();

    // Place the table view in the designer form
    this->ui->tableFrame->layout()->addWidget(this->m_tableView->GetWidget());

    // Mapper
    VTK_CREATE(vtkDataSetMapper, mapper);

    VTK_CREATE(vtkUnstructuredGrid, unstructuredGrid);
    mapper->SetInputData(unstructuredGrid);

    //mapper->ImmediateModeRenderingOn();

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
    this->ui->qvtkWidget->renderWindow()->AddRenderer(ren);

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
    this->m_connections->Connect(this->ui->qvtkWidget->renderWindow()->GetInteractor(),
        vtkCommand::KeyPressEvent, this, 
        SLOT(slotKeyPressed(vtkObject*, unsigned long, void*, void*)), 0, 1.0);

    m_mesher.setIntervals(999999, 999999);

    m_mesher.onUpdate = std::bind1st(std::mem_fun(&Mesh3D::update), this);

};

Mesh3D::~Mesh3D()
{
    // The smart pointers should clean up for up


}

void Mesh3D::slotKeyPressed(vtkObject *, unsigned long, void *, void *)
{

    vtkRenderWindowInteractor *rwi = this->ui->qvtkWidget->renderWindow()->GetInteractor();

    std::string key = rwi->GetKeySym();

    // Output the key that was pressed

    std::cout << "Pressed " << key << std::endl;

    if (key == "Right")
    {

        drawMesh();

        this->ui->qvtkWidget->renderWindow()->Render();

    } else if (key == "Left")
    {

        drawMesh();

        this->ui->qvtkWidget->renderWindow()->Render();

    } else if (key == "Escape")
    {

        m_mesher.stopMeshing();

        this->ui->qvtkWidget->renderWindow()->Render();

    } 
    else if (key == "c" || key == "C")
    {

        m_mesher.mesh().reset();
        drawMesh();

        this->ui->qvtkWidget->renderWindow()->Render();

    }
    else if (key == "m" || key == "M")
    {

        generateMesh(12, 12, 12);
        drawMesh();

        this->ui->qvtkWidget->renderWindow()->Render();

    } else if (key == "w" || key == "W")
    {

        vtkRenderer* ren = this->ui->qvtkWidget->renderWindow()->GetRenderers()->GetFirstRenderer();

        vtkActor* actor = ren->GetActors()->GetLastActor();

        actor->GetProperty()->SetRepresentationToWireframe();
        actor->GetProperty()->LightingOff();

        this->ui->qvtkWidget->renderWindow()->Render();

    }
    else if (key == "s" || key == "S")
    {

        vtkRenderer* ren = this->ui->qvtkWidget->renderWindow()->GetRenderers()->GetFirstRenderer();

        vtkActor* actor = ren->GetActors()->GetLastActor();

        actor->GetProperty()->SetRepresentationToSurface();
        actor->GetProperty()->LightingOn();

        this->ui->qvtkWidget->renderWindow()->Render();

    }
    else if (key == "a" || key == "A")
    {
        for (unsigned int i = 5; i < 10; ++i)
            for (unsigned int j = i; j < 10; ++j)
                for (unsigned int k = j; k < 10; ++k)
                {

                    unsigned int nu = i + 1;
                    unsigned int nv = j + 1;
                    unsigned int nw = k + 1;

                    bool res = generateMesh(nu, nv, nw);
                    drawMesh();

                    this->ui->qvtkWidget->renderWindow()->Render();

                    if (!res)
                        return;

                }
    }
    else if (key == "l" || key == "L")
    {

        bool res = loadMesh();
        drawMesh();

        this->ui->qvtkWidget->renderWindow()->Render();

    }

}

// Action to be taken upon file open 
void Mesh3D::slotOpenFile()
{

}

void Mesh3D::slotExit() 
{

    qApp->exit();

}

bool Mesh3D::generateMesh(const unsigned int nu, const unsigned int nv, const unsigned int nw)
{

    std::cout << "Start..." << std::endl;

    const double d = 0.25;

    std::cout << "nu = " << nu << std::endl;
    std::cout << "nv = " << nv << std::endl;
    std::cout << "nw = " << nw << std::endl;

    std::cout << "d = " << d << std::endl;

    CGeoCoordinate<double> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex2(nu * d, 0.0, 0.0);
    CGeoCoordinate<double> aVertex3(nu * d, nv * d, 0.0);
    CGeoCoordinate<double> aVertex4(0.0, nv * d, 0.0);
    CGeoCoordinate<double> aVertex5(0.0, 0.0, nw * d);
    CGeoCoordinate<double> aVertex6(nu * d, 0.0, nw * d);
    CGeoCoordinate<double> aVertex7(nu * d, nv * d, nw * d);
    CGeoCoordinate<double> aVertex8(0.0, nv * d, nw * d);

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

    aBasicMesher.generate(aHexahedron, nu, nv, nw, true);
    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CMshMesh<double> aBoundaryMesh = aMesh.extractBoundary(1E-6);

    aBoundaryMesh.generateFaces(1E-4);

    std::vector<CGeoCoordinate<double> > sInteriorPoints;

    QElapsedTimer timer;

    timer.start();

    bool res = false;

    try
    {

        // Volume mesh size should be larger than surface mesh size
        res = m_mesher.generate(aBoundaryMesh, 99999, sInteriorPoints, d, 0.1, 1E-12);

        std::cout << "Done. Elapsed time = " << timer.elapsed() / 1000.0 << " s" << std::endl;
        std::cout << "Number of elements: " << m_mesher.mesh().nbElements() << std::endl;

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
    aPosGmsh.save(aField, "tetra_volume.msh", "tetras");

    /*
    CMshMesh<double> aBoundaryMesh = m_mesher.mesh().extractBoundary(1E-4);
    aField.setMesh(aBoundaryMesh);
    aPosGmsh.save(aField, "tetra_surface.msh", "tris");
    */

    return res;

}

bool Mesh3D::loadMesh()
{

    std::cout << "Start..." << std::endl;

    std::cout << "fileName = " << m_fileName.toStdString() << std::endl;

    CPosGmsh<double> aPosGmsh;
    CPdeField<double> aField;

    aPosGmsh.load(aField, m_fileName.toStdString());

    aField.mesh().mergeNodes(0.01);
    aField.mesh().generateFaces(0.0001);

    CMshMesh<double> aBoundaryMesh = aField.mesh();

    std::vector<CGeoCoordinate<double> > sInteriorPoints;

    // Volume mesh size should be larger than surface mesh size
    bool res = m_mesher.generate(aBoundaryMesh, 200000, sInteriorPoints, 5, 0.1, 0.001);

    aField.setMesh(m_mesher.mesh());
    aPosGmsh.save(aField, "tetra_volume_load.msh", "tetras");

    std::cout << "Done." << std::endl;
    std::cout << "Number of elements: " << m_mesher.mesh().nbElements() << std::endl;

    return res;

}

void Mesh3D::drawMesh()
{

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

        if (anElement.elementType() == ET_TETRAHEDRON)
        {
            vtkIdType ptIds[] = {static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(0))),
                static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(1))),
                static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(2))),
                static_cast<vtkIdType> (m_mesher.mesh().nodeIndex(anElement.nodeId(3)))};

            unstructuredGrid->InsertNextCell(VTK_TETRA, 4, ptIds);
        }

    }

    vtkRenderer* ren = this->ui->qvtkWidget->renderWindow()->GetRenderers()->GetFirstRenderer();

    vtkMapper* mapper = ren->GetActors()->GetLastActor()->GetMapper();

    (reinterpret_cast<vtkDataSetMapper* >(mapper))->SetInputData(unstructuredGrid);

}

int Mesh3D::update(bool updateData)
{

    if (updateData)
    {

        m_mutex.lock();

        drawMesh();

        this->ui->qvtkWidget->renderWindow()->Render();

        m_mutex.unlock();

    }

    QCoreApplication::processEvents();

    return 0;

}

