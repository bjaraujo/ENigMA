// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "ui_TetrahedronIntersection.h"
#include "TetrahedronIntersection.h"

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
#include <vtkVertexGlyphFilter.h>

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// Constructor
TetrahedronIntersection::TetrahedronIntersection() 
{

    this->ui = new Ui_TetrahedronIntersection;
    this->ui->setupUi(this);

    // Qt Table View
    this->m_tableView = vtkSmartPointer<vtkQtTableView>::New();

    // Place the table view in the designer form
    this->ui->tableFrame->layout()->addWidget(this->m_tableView->GetWidget());

    // Mesh
    VTK_CREATE(vtkDataSetMapper, meshMapper);
    meshMapper->ImmediateModeRenderingOn();

    this->m_meshActor = vtkSmartPointer<vtkActor>::New();
    this->m_meshActor->SetMapper(meshMapper);

    m_meshActor->GetProperty()->SetRepresentationToSurface();
    m_meshActor->GetProperty()->LightingOn();

    // Point
    VTK_CREATE(vtkPolyDataMapper, pointMapper);
    pointMapper->ImmediateModeRenderingOn();

    this->m_pointActor = vtkSmartPointer<vtkActor>::New();
    this->m_pointActor->SetMapper(pointMapper);

    m_pointActor->GetProperty()->SetRepresentationToSurface();
    m_pointActor->GetProperty()->LightingOn();
    m_pointActor->GetProperty()->SetPointSize(5);

    // Tetrahedron
    VTK_CREATE(vtkDataSetMapper, tetrahedronMapper);
    tetrahedronMapper->ImmediateModeRenderingOn();

    this->m_tetrahedronActor = vtkSmartPointer<vtkActor>::New();
    this->m_tetrahedronActor->SetMapper(tetrahedronMapper);

    m_tetrahedronActor->GetProperty()->SetRepresentationToSurface();
    m_tetrahedronActor->GetProperty()->LightingOn();

    // Axes
    this->m_axesActor = vtkSmartPointer<vtkAxesActor>::New();

    // VTK Renderer
    VTK_CREATE(vtkRenderer, ren);

    // Add Actor to renderer
    ren->AddActor(m_axesActor);
    ren->AddActor(m_meshActor);
    ren->AddActor(m_pointActor);
    ren->AddActor(m_tetrahedronActor);

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

    // Create a PolyData
    buildMesh();

    drawMesh();
    drawPoint();
    drawTetrahedron();

};

TetrahedronIntersection::~TetrahedronIntersection()
{
    // The smart pointers should clean up for up


}

void TetrahedronIntersection::slotKeyPressed(vtkObject *, unsigned long, void *, void *)
{

    vtkRenderWindowInteractor *rwi = this->ui->qvtkWidget->GetRenderWindow()->GetInteractor();

    std::string key = rwi->GetKeySym();

    // Output the key that was pressed

    //std::cout << "Pressed " << key << std::endl;

    if (key == "Right")
    {

        m_point.z() += 0.01;

        drawPoint();
        drawTetrahedron();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    } else if (key == "Left")
    {

        m_point.z() -= 0.01;

        drawPoint();
        drawTetrahedron();

        this->ui->qvtkWidget->GetRenderWindow()->Render();

    } else if (key == "r" || key == "R")
    {

        drawMesh();
        drawPoint();

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
void TetrahedronIntersection::slotOpenFile()
{

}

void TetrahedronIntersection::slotExit() 
{

    qApp->exit();

}

void TetrahedronIntersection::buildMesh()
{

    std::vector<CMshNode<double> > sNodes;

    CMshNode<double> aNode1;
    CMshNode<double> aNode2;
    CMshNode<double> aNode3;
    CMshNode<double> aNode4;
    CMshNode<double> aNode5;
    CMshNode<double> aNode6;

    m_point << 0, 0, 0;

    aNode1 << 1.0, 0.875, 0.125;
    aNode2 << 0.875, 1.0, 0.125;
    aNode3 << 1.0, 1.0, 0.125;
    aNode4 << 1.0, 0.875, 0.25;
    aNode5 << 0.875, 1.0, 0.25;
    aNode6 << 1.0, 1.0, 0.25;

    sNodes.push_back(aNode1);
    sNodes.push_back(aNode2);
    sNodes.push_back(aNode3);
    sNodes.push_back(aNode4);
    sNodes.push_back(aNode5);
    sNodes.push_back(aNode6);

    for (unsigned int i = 0; i < sNodes.size(); ++i)
    {
        m_mesh.addNode(i + 1, sNodes[i]);

        m_point += sNodes[i];
    }

    m_point /= sNodes.size();

    CMshElement<double> anElement1(ET_TRIANGLE);
    CMshElement<double> anElement2(ET_TRIANGLE);
    CMshElement<double> anElement3(ET_TRIANGLE);
    CMshElement<double> anElement4(ET_TRIANGLE);
    CMshElement<double> anElement5(ET_TRIANGLE);
    CMshElement<double> anElement6(ET_TRIANGLE);
    CMshElement<double> anElement7(ET_TRIANGLE);
    CMshElement<double> anElement8(ET_TRIANGLE);

    // 1 - 2 5 3
    anElement1.addNodeId(2);
    anElement1.addNodeId(5);
    anElement1.addNodeId(3);
    m_mesh.addElement(1, anElement1);

    // 2 - 1 3 6
    anElement2.addNodeId(1);
    anElement2.addNodeId(3);
    anElement2.addNodeId(6);
    m_mesh.addElement(2, anElement2);

    // 3 - 3 5 6
    anElement3.addNodeId(3);
    anElement3.addNodeId(5);
    anElement3.addNodeId(6);
    m_mesh.addElement(3, anElement3);

    // 4 - 1 6 4
    anElement4.addNodeId(1);
    anElement4.addNodeId(6);
    anElement4.addNodeId(4);
    m_mesh.addElement(4, anElement4);

    // 5 - 6 5 4
    anElement5.addNodeId(6);
    anElement5.addNodeId(5);
    anElement5.addNodeId(4);
    m_mesh.addElement(5, anElement5);

    // 6 - 3 1 2
    anElement6.addNodeId(3);
    anElement6.addNodeId(1);
    anElement6.addNodeId(2);
    m_mesh.addElement(6, anElement6);

    // 7 - 4 5 2
    anElement7.addNodeId(4);
    anElement7.addNodeId(5);
    anElement7.addNodeId(2);
    m_mesh.addElement(7, anElement7);

    // 8 - 1 4 2
    anElement8.addNodeId(1);
    anElement8.addNodeId(4);
    anElement8.addNodeId(2);
    m_mesh.addElement(8, anElement8);

}

void TetrahedronIntersection::drawMesh()
{

    VTK_CREATE(vtkPoints, points);

    for (Integer i = 0; i < m_mesh.nbNodes(); ++i)
    {

        Integer aNodeId = m_mesh.nodeId(i);

        CMshNode<double> aNode = m_mesh.node(aNodeId);

        points->InsertNextPoint(aNode.x(), aNode.y(), aNode.z());

    }

    VTK_CREATE(vtkUnstructuredGrid, unstructuredGrid);

    unstructuredGrid->SetPoints(points);

    for (Integer i = 0; i < m_mesh.nbElements(); ++i)
    {

        Integer anElementId = m_mesh.elementId(i);

        CMshElement<double> anElement = m_mesh.element(anElementId);

        if (anElement.elementType() == ET_TRIANGLE)
        {            
            vtkIdType ptIds[] = {static_cast<vtkIdType> (m_mesh.nodeIndex(anElement.nodeId(0))),
                static_cast<vtkIdType> (m_mesh.nodeIndex(anElement.nodeId(1))),
                static_cast<vtkIdType> (m_mesh.nodeIndex(anElement.nodeId(2)))};

            unstructuredGrid->InsertNextCell(VTK_TRIANGLE, 3, ptIds);
        }

    }

    vtkMapper* mapper = m_meshActor->GetMapper();

    (reinterpret_cast<vtkDataSetMapper* >(mapper))->SetInputData(unstructuredGrid);

}

void TetrahedronIntersection::drawPoint()
{

    VTK_CREATE(vtkPoints, points);

    points->InsertNextPoint(m_point.x(), m_point.y(), m_point.z());

    VTK_CREATE(vtkPolyData, pointsPolydata);

    pointsPolydata->SetPoints(points);

    VTK_CREATE(vtkVertexGlyphFilter, vertexFilter);

    vertexFilter->SetInputData(pointsPolydata);

    vertexFilter->Update();

    // Setup colors
    unsigned char yellow[3] = {255, 255, 0};

    VTK_CREATE(vtkUnsignedCharArray, colors);

    colors->SetNumberOfComponents(3);
    colors->SetName ("Colors");
    colors->InsertNextTypedTuple(yellow);

    VTK_CREATE(vtkPolyData, polydata);

    polydata->ShallowCopy(vertexFilter->GetOutput());

    polydata->GetPointData()->SetScalars(colors);

    vtkMapper* mapper = m_pointActor->GetMapper();

    (reinterpret_cast<vtkPolyDataMapper* >(mapper))->SetInputData(polydata);

}

void TetrahedronIntersection::drawTetrahedron()
{

    CGeoTetrahedron<double> aTetrahedron;

    for (unsigned int i = 0; i < 3; ++i)
    {

        unsigned int aNodeId = m_mesh.nodeId(i);
        CMshNode<double> aNode = m_mesh.node(aNodeId);

        CGeoCoordinate<double> aVertex = aNode;
        aTetrahedron.addVertex(aVertex);

    }

    CGeoCoordinate<double> aVertex = m_point;
    aTetrahedron.addVertex(aVertex);

    // Setup colors
    unsigned char activeColor[3];
    unsigned char red[3] = {255, 0, 0};
    unsigned char green[3] = {0, 255, 0};

    activeColor[0] = green[0];
    activeColor[1] = green[1];
    activeColor[2] = green[2];

    bool bOK = true;
    double tol = 1E-5;

    for (Integer i = 3; i < m_mesh.nbNodes(); ++i)
    {

        unsigned int aNodeId = m_mesh.nodeId(i);

        CMshNode<double> aNode = m_mesh.node(aNodeId);

        CGeoCoordinate<double> aPoint = aNode;

        if (aTetrahedron.contains(aPoint, tol))
        {
            bOK = false;
            break;
        }

    }

    CGeoTriangle<double> aTriangle1;

    aTriangle1.addVertex(aTetrahedron.vertex(0));
    aTriangle1.addVertex(aTetrahedron.vertex(1));
    aTriangle1.addVertex(aTetrahedron.vertex(3));

    CGeoTriangle<double> aTriangle2;

    aTriangle2.addVertex(aTetrahedron.vertex(1));
    aTriangle2.addVertex(aTetrahedron.vertex(2));
    aTriangle2.addVertex(aTetrahedron.vertex(3));

    CGeoTriangle<double> aTriangle3;

    aTriangle3.addVertex(aTetrahedron.vertex(2));
    aTriangle3.addVertex(aTetrahedron.vertex(0));
    aTriangle3.addVertex(aTetrahedron.vertex(3));

    for (Integer i = 3; i < m_mesh.nbElements(); ++i)
    {

        Integer anElementId = m_mesh.elementId(i);

        CMshElement<double> anElement = m_mesh.element(anElementId);

        if (anElement.elementType() == ET_TRIANGLE)
        {            

            CGeoTriangle<double> aTriangle;

            for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
            {

                CGeoCoordinate<double> aVertex;

                aVertex = m_mesh.node(anElement.nodeId(j));

                aTriangle.addVertex(aVertex);

            }

            CGeoIntersectionType anIntType;

            if (aTriangle.intersects(aTriangle1, anIntType, tol))
            {
                if (anIntType == IT_INTERNAL)
                {
                    bOK = false;
                    break;
                }
            }

            if (aTriangle.intersects(aTriangle2, anIntType, tol))
            {
                if (anIntType == IT_INTERNAL)
                {
                    bOK = false;
                    break;
                }
            }

            if (aTriangle.intersects(aTriangle3, anIntType, tol))
            {
                if (anIntType == IT_INTERNAL)
                {
                    bOK = false;
                    break;
                }
            }

        }

    }

    if (!bOK)
    {
        activeColor[0] = red[0];
        activeColor[1] = red[1];
        activeColor[2] = red[2];
    }

    if (bOK)
        std::cout << "Tetra: OK!" << std::endl;
    else
        std::cout << "Tetra: Not OK!" << std::endl;

    VTK_CREATE(vtkUnsignedCharArray, colors);

    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    VTK_CREATE(vtkPoints, points);

    for (unsigned int i = 0; i < 4; ++i)
    {

        CGeoCoordinate<double> aVertex = aTetrahedron.vertex(i);
        points->InsertNextPoint(aVertex.x(), aVertex.y(), aVertex.z());

        colors->InsertNextTypedTuple(activeColor);

    }

    VTK_CREATE(vtkUnstructuredGrid, unstructuredGrid);

    unstructuredGrid->SetPoints(points);

    vtkIdType ptIds[] = {0, 1, 2, 3};

    unstructuredGrid->InsertNextCell(VTK_TETRA, 4, ptIds);

    unstructuredGrid->GetPointData()->SetScalars(colors);

    vtkMapper* mapper = m_tetrahedronActor->GetMapper();

    (reinterpret_cast<vtkDataSetMapper* >(mapper))->SetInputData(unstructuredGrid);

}
