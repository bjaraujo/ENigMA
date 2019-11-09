// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "ui_TriangleIntersection.h"
#include "TriangleIntersection.h"

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
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// Constructor
TriangleIntersection::TriangleIntersection() 
{

    this->ui = new Ui_TriangleIntersection;
    this->ui->setupUi(this);

    // Qt Table View
    this->m_tableView = vtkSmartPointer<vtkQtTableView>::New();

    // Place the table view in the designer form
    this->ui->tableFrame->layout()->addWidget(this->m_tableView->GetWidget());

    // Mapper
    VTK_CREATE(vtkPolyDataMapper, mapper);
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

    m_position << 0.0, 0.0, 0.0;

    m_p1.x() = (rand() % 100) / 100.0;
    m_p1.y() = (rand() % 100) / 100.0;
    m_p1.z() = (rand() % 100) / 100.0;

    m_p2.x() = (rand() % 100) / 100.0;
    m_p2.y() = (rand() % 100) / 100.0;
    m_p2.z() = (rand() % 100) / 100.0;

    m_p3.x() = (rand() % 100) / 100.0;
    m_p3.y() = (rand() % 100) / 100.0;
    m_p3.z() = (rand() % 100) / 100.0;

    m_p4.x() = (rand() % 100) / 100.0;
    m_p4.y() = (rand() % 100) / 100.0;
    m_p4.z() = (rand() % 100) / 100.0;

    m_p5.x() = (rand() % 100) / 100.0;
    m_p5.y() = (rand() % 100) / 100.0;
    m_p5.z() = (rand() % 100) / 100.0;

    m_p6.x() = (rand() % 100) / 100.0;
    m_p6.y() = (rand() % 100) / 100.0;
    m_p6.z() = (rand() % 100) / 100.0;

    // Create a PolyData
    drawTriangles();

};

TriangleIntersection::~TriangleIntersection()
{
    // The smart pointers should clean up for up


}

void TriangleIntersection::slotKeyPressed(vtkObject *, unsigned long, void *, void *)
{

    vtkRenderWindowInteractor *rwi = this->ui->qvtkWidget->GetRenderWindow()->GetInteractor();
    
    std::string key = rwi->GetKeySym();
 
    // Output the key that was pressed
     
    //std::cout << "Pressed " << key << std::endl;

    if (key == "Up")
    {

        m_position.z() += 0.01;

        drawTriangles();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    } else if (key == "Down")
    {

        m_position.z() -= 0.01;

        drawTriangles();

        this->ui->qvtkWidget->GetRenderWindow()->Render();
    }
    if (key == "Right")
    {

        m_position.x() += 0.01;

        drawTriangles();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    } else if (key == "Left")
    {

        m_position.x() -= 0.01;

        drawTriangles();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    } else if (key == "r" || key == "R")
    {

        m_p1.x() = (rand() % 100) / 100.0;
        m_p1.y() = (rand() % 100) / 100.0;
        m_p1.z() = (rand() % 100) / 100.0;

        m_p2.x() = (rand() % 100) / 100.0;
        m_p2.y() = (rand() % 100) / 100.0;
        m_p2.z() = (rand() % 100) / 100.0;

        m_p3.x() = (rand() % 100) / 100.0;
        m_p3.y() = (rand() % 100) / 100.0;
        m_p3.z() = (rand() % 100) / 100.0;

        m_p4.x() = (rand() % 100) / 100.0;
        m_p4.y() = (rand() % 100) / 100.0;
        m_p4.z() = (rand() % 100) / 100.0;

        m_p5.x() = (rand() % 100) / 100.0;
        m_p5.y() = (rand() % 100) / 100.0;
        m_p5.z() = (rand() % 100) / 100.0;

        m_p6.x() = (rand() % 100) / 100.0;
        m_p6.y() = (rand() % 100) / 100.0;
        m_p6.z() = (rand() % 100) / 100.0;

        drawTriangles();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    } else if (key == "w" || key == "W")
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
void TriangleIntersection::slotOpenFile()
{

}

void TriangleIntersection::slotExit() 
{
    
    qApp->exit();

}

void TriangleIntersection::drawTriangles()
{
    
    std::vector<CGeoCoordinate<double> > sVertices;

    CGeoCoordinate<double> aVertex1(+7.375, -12.800, -8.499);
    CGeoCoordinate<double> aVertex2(+8.225, -12.670, -7.275);
    CGeoCoordinate<double> aVertex3(+8.650, -13.730, -8.500);

    CGeoCoordinate<double> aVertex4(+8.649, -12.96, -8.500);
    CGeoCoordinate<double> aVertex5(+8.650, -13.73, -8.500);
    CGeoCoordinate<double> aVertex6(+7.812, -13.68, -8.500);

    /*
    CGeoCoordinate<double> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex2(1.0, 1.0, 0.0);
    CGeoCoordinate<double> aVertex3(1.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex4(0.0 + m_position.z(), 0.0, 0.0);
    CGeoCoordinate<double> aVertex5(0.0 + m_position.z(), -1.0, 0.0);
    CGeoCoordinate<double> aVertex6(1.0 + m_position.z(), 0.0, 0.0);
    */

    sVertices.push_back(aVertex1);
    sVertices.push_back(aVertex2);
    sVertices.push_back(aVertex3);
    sVertices.push_back(aVertex4);
    sVertices.push_back(aVertex5);
    sVertices.push_back(aVertex6);

    m_triangle1.reset();
    m_triangle1.addVertex(sVertices[0]);
    m_triangle1.addVertex(sVertices[1]);
    m_triangle1.addVertex(sVertices[2]);

    m_triangle2.reset();
    m_triangle2.addVertex(sVertices[3]);
    m_triangle2.addVertex(sVertices[4]);
    m_triangle2.addVertex(sVertices[5]);

    // Setup colors
    unsigned char red[3] = {255, 0, 0};
    unsigned char yellow[3] = {255, 255, 0};
    unsigned char cyan[3] = {0, 255, 255};
 
    VTK_CREATE(vtkUnsignedCharArray, colors);

    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    VTK_CREATE(vtkPoints, points);

    CGeoIntersectionType anIntersectionType;

    if(m_triangle1.intersects(m_triangle2, anIntersectionType, 1E-2))
        std::cout << "-> INTERSECTS" << std::endl;
    else
        std::cout << "-> DOES NOT INTERSECT" << std::endl;

    if (anIntersectionType == IT_VERTEX)
        std::cout << "--> IT_VERTEX" << std::endl;
    else if (anIntersectionType == IT_EDGE)
        std::cout << "--> IT_EDGE" << std::endl;
    else if (anIntersectionType == IT_COINCIDENT)
        std::cout << "--> IT_COINCIDENT" << std::endl;
    else if (anIntersectionType == IT_INTERNAL)
        std::cout << "--> IT_INTERNAL" << std::endl;
    else if (anIntersectionType == IT_SWAP)
        std::cout << "--> IT_SWAP" << std::endl;
    else if (anIntersectionType == IT_NONE)
        std::cout << "--> IT_NONE" << std::endl;

    for (unsigned int i = 0; i < sVertices.size(); ++i)
    {
        points->InsertNextPoint(sVertices[i].x(), sVertices[i].y(), sVertices[i].z());

        if (i <= 2)
        {

            if (anIntersectionType == IT_NONE)
                colors->InsertNextTypedTuple(cyan);
            else
                colors->InsertNextTypedTuple(red);

        }
        else
        {
            
            if (anIntersectionType == IT_NONE)
                colors->InsertNextTypedTuple(yellow);
            else
                colors->InsertNextTypedTuple(red);

        }

    }

    VTK_CREATE(vtkTriangle, triangle1);

    triangle1->GetPointIds()->SetId(0, 0);
    triangle1->GetPointIds()->SetId(1, 1);
    triangle1->GetPointIds()->SetId(2, 2);

    VTK_CREATE(vtkTriangle, triangle2);

    triangle2->GetPointIds()->SetId(0, 3);
    triangle2->GetPointIds()->SetId(1, 4);
    triangle2->GetPointIds()->SetId(2, 5);

    VTK_CREATE(vtkCellArray, triangles);
    triangles->InsertNextCell(triangle1);
    triangles->InsertNextCell(triangle2);

    VTK_CREATE(vtkPolyData, polygonPolyData);

    polygonPolyData->SetPoints(points);
    polygonPolyData->SetPolys(triangles);

    polygonPolyData->GetPointData()->SetScalars(colors);

    vtkRenderer* ren = this->ui->qvtkWidget->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

    vtkMapper* mapper = ren->GetActors()->GetLastActor()->GetMapper();

    (reinterpret_cast<vtkPolyDataMapper* >(mapper))->SetInputData(polygonPolyData);

}
