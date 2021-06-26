// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "ui_ClipPolyhedron.h"
#include "ClipPolyhedron.h"

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
#include <vtkUnsignedCharArray.h>

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// Constructor
ClipPolyhedron::ClipPolyhedron() 
{
    this->ui = new Ui_ClipPolyhedron;
    this->ui->setupUi(this);

    this->setWindowTitle("Polyhedron");

    // Qt Table View
    this->m_tableView = vtkSmartPointer<vtkQtTableView>::New();

    // Place the table view in the designer form
    this->ui->tableFrame->layout()->addWidget(this->m_tableView->GetWidget());

    // Mapper
    VTK_CREATE(vtkPolyDataMapper, mapper);
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
    connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));

    m_connections = vtkEventQtSlotConnect::New();
    this->m_connections->Connect(this->ui->qvtkWidget->renderWindow()->GetInteractor(),
        vtkCommand::KeyPressEvent, this, 
        SLOT(slotKeyPressed(vtkObject*, unsigned long, void*, void*)), 0, 1.0);

    // Plane equation
    m_plane.set(CGeoNormal<double>(1.0, 0.0, 0.0), 0.5);

    m_position << 0.0, 0.0, 0.0;

    std::cout << "Cutting plane normal: " << m_plane.normal() << std::endl;

    m_fraction = 0.5;

    m_size = 1.0;

    m_errors = 0;

    // Create a PolyData
    drawClippedCube();

};

ClipPolyhedron::~ClipPolyhedron()
{
    // The smart pointers should clean up for up


}

void ClipPolyhedron::slotKeyPressed(vtkObject *, unsigned long, void *, void *)
{

    vtkRenderWindowInteractor *rwi = this->ui->qvtkWidget->renderWindow()->GetInteractor();

    std::string key = rwi->GetKeySym();

    // Output the key that was pressed

    //std::cout << "Pressed " << key << std::endl;

    if (key == "Right")
    {

        m_fraction += 0.01;
        m_fraction = std::min(m_fraction, 1.0);

        m_plane.setD(m_plane.d() + 0.01);

        std::cout << "d = " << m_plane.d() << std::endl;

        drawClippedCube();

        this->ui->qvtkWidget->interactor()->Render();

    } else if (key == "Left")
    {

        m_fraction -= 0.01;
        m_fraction = std::max(m_fraction, 0.0);

        m_plane.setD(m_plane.d() - 0.01);
        
        std::cout << "d = " << m_plane.d() << std::endl;

        drawClippedCube();

        this->ui->qvtkWidget->interactor()->Render();

    } else if (key == "r" || key == "R")
    {

        m_plane.normal().x() = (rand() % 200) / 100.0 - 1.0;
        m_plane.normal().y() = (rand() % 200) / 100.0 - 1.0;
        m_plane.normal().z() = (rand() % 200) / 100.0 - 1.0;

        std::cout << "Cutting plane normal: " << m_plane.normal() << std::endl;

        drawClippedCube();

        this->ui->qvtkWidget->interactor()->Render();

    } else if (key == "w" || key == "W")
    {

        vtkRenderer* ren = this->ui->qvtkWidget->renderWindow()->GetRenderers()->GetFirstRenderer();

        vtkActor* actor = ren->GetActors()->GetLastActor();

        actor->GetProperty()->SetRepresentationToWireframe();
        actor->GetProperty()->LightingOff();

        this->ui->qvtkWidget->interactor()->Render();

    }
    else if (key == "s" || key == "S")
    {

        vtkRenderer* ren = this->ui->qvtkWidget->renderWindow()->GetRenderers()->GetFirstRenderer();

        vtkActor* actor = ren->GetActors()->GetLastActor();

        actor->GetProperty()->SetRepresentationToSurface();
        actor->GetProperty()->LightingOn();

        this->ui->qvtkWidget->interactor()->Render();

    }
}

void ClipPolyhedron::slotExit() 
{

    qApp->exit();

}

void ClipPolyhedron::drawClippedCube()
{

    // Build polyhedron
    std::vector<CGeoCoordinate<double> > sVertices;

    sVertices.push_back(CGeoCoordinate<double>(0, 0, m_size));
    sVertices.push_back(CGeoCoordinate<double>(m_size, 0, m_size));
    sVertices.push_back(CGeoCoordinate<double>(m_size, m_size, m_size));
    sVertices.push_back(CGeoCoordinate<double>(0, m_size, m_size));
    sVertices.push_back(CGeoCoordinate<double>(0, 0, 0));
    sVertices.push_back(CGeoCoordinate<double>(m_size, 0, 0));
    sVertices.push_back(CGeoCoordinate<double>(m_size, m_size, 0));
    sVertices.push_back(CGeoCoordinate<double>(0, m_size, 0));

    CGeoHexahedron<double> aHexahedron;

    for (Integer i = 0; i < 8; ++i)
    {
        CGeoCoordinate<double> aVertex;

        aVertex = m_position + sVertices[i];
        aHexahedron.addVertex(aVertex);

    }

    CGeoPolyhedron<double> aPolyhedron(aHexahedron);

    Integer iter;
    double frac_req, frac_act;

    frac_req = m_fraction;

    CGeoPolyhedron<double> aNewPolyhedron;
    CGeoPolygon<double> aNewPolygon;
    Integer aNewPolygonId = 999;

    aPolyhedron = aPolyhedron.clip(aNewPolygon, aNewPolygonId, m_plane, frac_req, frac_act, iter, 100, 1E-5, 1E-5);
    //aPolyhedron = aPolyhedron.cut(aNewPolyhedron, aNewPolygon, aNewPolygonId, m_plane, frac_req, frac_act, iter, 100, 1E-5, 1E-5);

    if (fabs(frac_req - frac_act) > 1E-6)
        m_errors++;

    std::cout << "************* Volume fraction *************" << std::endl;
    std::cout << "Required: " << frac_req << std::endl;
    std::cout << "Actual: " << frac_act << std::endl;
    std::cout << "Iterations: " << iter << std::endl;
    std::cout << "Error count: " << m_errors << std::endl;

    aPolyhedron.calculateSurfaceArea();
    aPolyhedron.calculateVolume();

    // Draw clipped polyhedron
    VTK_CREATE(vtkPoints, points);
    VTK_CREATE(vtkCellArray, polygons);

    // Setup colors
    unsigned char yellow[3] = {255, 255, 0};
    unsigned char cyan[3] = {0, 255, 255};

    VTK_CREATE(vtkUnsignedCharArray, colors);

    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    int np = 0;

    // Clipped polygon 

    for (Integer i = 0; i < aPolyhedron.nbPolygons(); ++i)
    {
        VTK_CREATE(vtkPolygon, polygon);

        Integer id = aPolyhedron.polygonId(i);

        polygon->GetPointIds()->SetNumberOfIds(aPolyhedron.polygon(id).polyline().nbVertices() - 1);

        for (Integer j = 0; j < aPolyhedron.polygon(id).polyline().nbVertices() - 1; ++j)
        {

            CGeoCoordinate<double> vertex = aPolyhedron.polygon(id).polyline().vertex(j);

            polygon->GetPointIds()->SetId(j, np);

            points->InsertNextPoint(vertex.x(), vertex.y(), vertex.z());

            colors->InsertNextTypedTuple(yellow);

            np++;

        }

        // Add the polygon to a list of polygons
        polygons->InsertNextCell(polygon);

    }

    std::cout << ">>>>>>>>>>>>>>>>> " << aNewPolygon.polyline().nbLines() << std::endl;

    {
        VTK_CREATE(vtkPolygon, polygon);

        polygon->GetPointIds()->SetNumberOfIds(aNewPolygon.polyline().nbVertices() - 1);

        for (Integer j = 0; j < aNewPolygon.polyline().nbVertices() - 1; ++j)
        {
            CGeoCoordinate<double> vertex = aNewPolygon.polyline().vertex(j);

            polygon->GetPointIds()->SetId(j, np);

            points->InsertNextPoint(vertex.x(), vertex.y(), vertex.z());

            colors->InsertNextTypedTuple(cyan);

            np++;

        }

        // Add the polygon to a list of polygons
        polygons->InsertNextCell(polygon);
    }

    // New polygon 
    /*
    bool bShowNewPolygon = true;
    if (bShowNewPolygon)
    {

        for (Integer i = 0; i < aNewPolyhedron.nbPolygons(); ++i)
        {

            VTK_CREATE(vtkPolygon, polygon);

            Integer id = aNewPolyhedron.polygonId(i);

            polygon->GetPointIds()->SetNumberOfIds(aNewPolyhedron.polygon(id).polyline().nbVertices() - 1);

            for (Integer j = 0; j < aNewPolyhedron.polygon(id).polyline().nbVertices() - 1; ++j)
            {

                CGeoCoordinate<double> vertex = aNewPolyhedron.polygon(id).polyline().vertex(j);

                polygon->GetPointIds()->SetId(j, np);

                points->InsertNextPoint(vertex.x(), vertex.y(), vertex.z());

                colors->InsertNextTypedTuple(cyan);

                np++;

            }

            // Add the polygon to a list of polygons
            polygons->InsertNextCell(polygon);

        }

    }
    */

    VTK_CREATE(vtkPolyData, polygonPolyData);

    polygonPolyData->SetPoints(points);
    polygonPolyData->SetPolys(polygons);

    polygonPolyData->GetPointData()->SetScalars(colors);

    vtkRenderer* ren = this->ui->qvtkWidget->renderWindow()->GetRenderers()->GetFirstRenderer();

    vtkMapper* mapper = ren->GetActors()->GetLastActor()->GetMapper();

    (reinterpret_cast<vtkPolyDataMapper* >(mapper))->SetInputData(polygonPolyData);
}
