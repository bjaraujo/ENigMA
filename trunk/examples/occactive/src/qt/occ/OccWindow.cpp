
#include <QtGui/QInputDialog>
#include <QtCore/QDebug>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp_Explorer.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Curve.hxx>
#include <Geom2d_Curve.hxx>

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ShapeAnalysis_Curve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepTools.hxx>
#include <BRepAdaptor_Curve.hxx>

#include <BRepOffsetAPI_Sewing.hxx>
#include <Handle_ShapeBuild_ReShape.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <Prs3d_ShapeTool.hxx>

#include <Geom_CartesianPoint.hxx>
#include <GeomLProp_CLProps.hxx>
#include <GeomLProp_SLProps.hxx>

#include <TShort_HArray1OfShortReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <Poly_Array1OfTriangle.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>

#include <TColStd_ListOfInteger.hxx>

#include <AIS_Point.hxx>
#include <AIS_Line.hxx>
#include <AIS_Triangulation.hxx>
#include <AIS_Shape.hxx>

#include <Prs3d_Drawer.hxx>
#include <Prs3d_LineAspect.hxx>
#include <Prs3d_ShadingAspect.hxx>

#include <Graphic3d_AspectFillArea3d.hxx>

#include "OccWindow.h"

OccWindow::OccWindow()
{

    m_displayConnection = new Aspect_DisplayConnection();
    m_graphicDriver = new OpenGl_GraphicDriver(m_displayConnection);

    TCollection_ExtendedString aName("3DV");

    m_viewer = new V3d_Viewer(m_graphicDriver, aName.ToExtString());

    m_viewer->CreateView();
    m_viewer->SetDefaultLights();
    m_viewer->SetLightOn();

    m_viewer->PrivilegedPlane().SetDirection(gp_Dir(0.0, 0.0, 1.0));
    m_viewer->PrivilegedPlane().SetLocation(gp_Pnt(0.0, 0.0, 0.0));

    m_context = new AIS_InteractiveContext(m_viewer);

    m_context->SetDisplayMode(AIS_WireFrame);

    m_context->SetHilightColor(Quantity_NOC_CYAN2);
    m_context->SelectionColor(Quantity_NOC_CYAN4);
    m_context->SetPreselectionColor(Quantity_NOC_CYAN4);

    m_occWidget.setParent(this);
    m_occWidget.setContext(m_context);

    this->setWidget(&m_occWidget);

    this->setMinimumWidth(640);
    this->setMinimumHeight(480);

    connect(&m_occWidget, SIGNAL(onSelectionChanged()), this, SLOT(selectionChanged()));

}

OccWindow::~OccWindow()
{

    // Release memory 
    // According to Open CASCADE Technology 6.7.0 - Coding Rules, page 14
    m_displayConnection.Nullify();
    m_graphicDriver.Nullify();
    m_viewer.Nullify();
    m_context.Nullify();

}

void OccWindow::selectionChanged()
{

    emit onSelectionChanged();

}

void OccWindow::eraseAll()
{

    try
    {

        m_context->CloseAllContexts();
        m_context->RemoveAll();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void OccWindow::selectVertices()
{

    // Select vertices
    m_context->OpenLocalContext(); 
    m_context->ActivateStandardMode(TopAbs_VERTEX); 

}

void OccWindow::selectEdges()
{

    // Select edges
    m_context->OpenLocalContext(); 
    m_context->ActivateStandardMode(TopAbs_EDGE); 

}

void OccWindow::selectWires()
{

    // Select edges
    m_context->OpenLocalContext(); 
    m_context->ActivateStandardMode(TopAbs_WIRE); 

}

void OccWindow::selectFaces()
{

    // Select faces
    m_context->OpenLocalContext(); 
    m_context->ActivateStandardMode(TopAbs_FACE); 

}

void OccWindow::selectSolids()
{

    // Select solids
    m_context->OpenLocalContext(); 
    m_context->ActivateStandardMode(TopAbs_SOLID); 

}

void OccWindow::selectAll()
{

    // Select all
    m_context->CloseAllContexts();

}

void OccWindow::viewWireframe()
{

    for(m_context->InitCurrent(); m_context->MoreCurrent(); m_context->NextCurrent())
        m_context->SetDisplayMode(m_context->Current(), AIS_WireFrame);

}

void OccWindow::viewShaded()
{

    for(m_context->InitCurrent(); m_context->MoreCurrent(); m_context->NextCurrent())
        m_context->SetDisplayMode(m_context->Current(), AIS_Shaded);
    
}

void OccWindow::viewFront()
{
    
    m_occWidget.front();
    m_occWidget.fitAll();

}

void OccWindow::viewBack()
{
    
    m_occWidget.back();
    m_occWidget.fitAll();

}

void OccWindow::viewLeft()
{
    
    m_occWidget.left();
    m_occWidget.fitAll();

}

void OccWindow::viewRight()
{
    
    m_occWidget.right();
    m_occWidget.fitAll();

}

void OccWindow::viewTop()
{
    
    m_occWidget.top();
    m_occWidget.fitAll();

}

void OccWindow::viewBottom()
{
    
    m_occWidget.bottom();
    m_occWidget.fitAll();

}

void OccWindow::viewIsometric()
{
    
    m_occWidget.axo();
    m_occWidget.fitAll();

}    

void OccWindow::setColor(QColor& aColor)
{
            
    for (m_context->InitCurrent(); m_context->MoreCurrent(); m_context->NextCurrent())
        m_context->SetColor(m_context->Current(), Quantity_Color(aColor.red()/255., aColor.green()/255., aColor.blue()/255., Quantity_TOC_RGB));

}

void OccWindow::setMaterial(eMaterial aMaterial)
{

    Graphic3d_NameOfMaterial cMaterial;

    if (aMaterial == MAT_BRASS)
        cMaterial = Graphic3d_NOM_BRASS;
    else if (aMaterial == MAT_BRONZE)
        cMaterial = Graphic3d_NOM_BRONZE;
    else if (aMaterial == MAT_COPPER)
        cMaterial = Graphic3d_NOM_COPPER;
    else if (aMaterial == MAT_GOLD)
        cMaterial = Graphic3d_NOM_GOLD;
    else if (aMaterial == MAT_PEWTER)
        cMaterial = Graphic3d_NOM_PEWTER;
    else if (aMaterial == MAT_PLASTIC)
        cMaterial = Graphic3d_NOM_PLASTIC;
    else if (aMaterial == MAT_SILVER)
        cMaterial = Graphic3d_NOM_SILVER;
    else if (aMaterial == MAT_STEEL)
        cMaterial = Graphic3d_NOM_STEEL;

    for (m_context->InitCurrent(); m_context->MoreCurrent(); m_context->NextCurrent())
        m_context->SetMaterial(m_context->Current(), cMaterial);

}

void OccWindow::fitAll()
{

    m_occWidget.fitAll();

}

Handle(AIS_InteractiveContext) OccWindow::context()
{

    return m_context;

}
