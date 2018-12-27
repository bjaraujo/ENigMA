
#include <QtCore/QCoreApplication>
#include <QtCore/QFileInfo>
#include <QtCore/QDebug>

#include <AIS_InteractiveContext.hxx>
#include <AIS_InteractiveObject.hxx>
#include <AIS_Shape.hxx>
#include <AIS_Line.hxx>
#include <AIS_Triangulation.hxx>

#include <Graphic3d_AspectFillArea3d.hxx>

#include <Poly_Triangulation.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Builder.hxx>

#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepFill_Filling.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <GeomConvert.hxx>
#include <Geom_OffsetSurface.hxx>

#include <Poly_Triangle.hxx>

#include <Prs3d_LineAspect.hxx>
#include <Prs3d_ShadingAspect.hxx>
#include <Prs3d_ShapeTool.hxx>

#include <BRepClass3d.hxx>

#include <ShapeAnalysis_FreeBounds.hxx>
#include <ShapeAnalysis_Shell.hxx>
#include <ShapeFix_ShapeTolerance.hxx>

#include "OccManager.h"

OccManager::OccManager() : m_selectionModel(&m_treeViewModel)
{

    connect(&m_occWindow, SIGNAL(onSelectionChanged()), this, SLOT(selectionChanged()));

    connect(&m_selectionModel, SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this, SLOT(treeViewSelectionChanged(const QItemSelection&, const QItemSelection&)));

    this->reset();

}

OccManager::~OccManager()
{



}


OccWindow* OccManager::window() 
{ 

    return &m_occWindow; 

}

QStandardItemModel* OccManager::treeViewModel()
{

    return &m_treeViewModel;

}

QItemSelectionModel* OccManager::treeViewSelectionModel()
{

    return &m_selectionModel;

}

void OccManager::selectionChanged()
{

    try
    {

        emit onSelectionChanged();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << tr("Unknown exception:") << QString(this->metaObject()->className());
    }

}

void OccManager::reset()
{

    try
    {

        m_occWindow.eraseAll();

        m_treeViewModel.clear();

        m_headerItem = new QStandardItem();
        m_headerItem->setText(tr("Structure"));

        m_treeViewModel.setHorizontalHeaderItem(0, m_headerItem);

        m_modelRoot = new QStandardItem();
        m_modelRoot->setText(tr("Model"));
        m_treeViewModel.appendRow(m_modelRoot);

        m_mesher.reset();

        m_shapeMap.clear();

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

void OccManager::displayShape(const TopoDS_Shape& aShape)
{

    try
    {

        Handle(AIS_Shape) aNewShape = new AIS_Shape(aShape);

        aNewShape->Attributes()->SetFaceBoundaryDraw(Standard_True);
        aNewShape->Attributes()->FaceBoundaryAspect()->SetColor(Quantity_NOC_BLACK);
        aNewShape->Attributes()->FaceBoundaryAspect()->SetTypeOfLine(Aspect_TOL_SOLID);
        aNewShape->Attributes()->FaceBoundaryAspect()->SetWidth(2.0);

        m_occWindow.context()->Display(aNewShape, false);

        m_occWindow.context()->SetMaterial(aNewShape, Graphic3d_NOM_PLASTIC);
        m_occWindow.context()->SetColor(aNewShape, Quantity_Color(0.9, 0.9, 0.75, Quantity_TOC_RGB));

        m_occWindow.context()->SetDisplayMode(aNewShape, AIS_Shaded);

        aNewShape.Nullify();

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

void OccManager::importGeometry(int nFormat, QString& aFileName)
{

    try
    {

        m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeOff();
        m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetInteriorStyle(Aspect_IS_SOLID);
        m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeLineType(Aspect_TOL_SOLID);

        Handle(TopTools_HSequenceOfShape) sShapes = m_translate.importModel(nFormat, aFileName.toUtf8().constData());

        if (sShapes.IsNull() || sShapes->Length() == 0)
            return;

        for (int i = 1; i <= sShapes->Length(); ++i)
        {

            TopoDS_Shape aShape = sShapes->Value(i);

            this->displayShape(aShape);

            QStandardItem* aShapeItem = new QStandardItem(QString("Shape %0").arg(i - 1));

            m_modelRoot->appendRow(aShapeItem);

            m_shapeMap[aShapeItem->index()] = aShape;

            Integer j = 0;

            for (TopExp_Explorer aSurfExp(sShapes->Value(i), TopAbs_FACE); aSurfExp.More(); aSurfExp.Next())
            {

                TopoDS_Face aFace = TopoDS::Face(aSurfExp.Current());

                QStandardItem* faceItem = new QStandardItem(QString("Face %0").arg(++j));

                aShapeItem->appendRow(faceItem);

                m_shapeMap[faceItem->index()] = aFace;

            }

        }

        m_occWindow.context()->UpdateCurrentViewer();

        // Fit all
        m_occWindow.fitAll();

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

void OccManager::exportGeoemtry(int nFormat, QString& aFileName)
{

}

void OccManager::importMesh(int nFormat, QString& aFileName)
{

    m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeOn();
    m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetInteriorStyle(Aspect_IS_EMPTY);
    m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeColor(Quantity_NOC_DARKGREEN);
    m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeLineType(Aspect_TOL_SOLID);
    m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeWidth(1.0);

    if (nFormat == OccTranslator::FormatGmsh)
    {

        CPdeField<double> aField;
        CPosGmsh<double> aPosGmsh;

        aPosGmsh.load(aField, aFileName.toStdString());

        m_mesher.mesh().reset();
        m_mesher.mesh().addMesh(aField.mesh());
        m_mesher.mesh().invert();

    }
    else if(nFormat == OccTranslator::FormatSTL)
    {

        CStlUtils<double> aStlUtils;

        aStlUtils.load(aFileName.toStdString());

        m_mesher.mesh().reset();
        m_mesher.mesh().addMesh(aStlUtils.mesh());
        m_mesher.mesh().invert();

    }
    else
        return;

    std::vector<gp_Pnt> pnts;
    std::vector<Poly_Triangle> tris;

    for (Integer i = 0; i < m_mesher.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = m_mesher.mesh().nodeId(i);

        CMshNode<double>& aNode = m_mesher.mesh().node(aNodeId);

        gp_Pnt aPoint(aNode.x(), aNode.y(),aNode.z());

        pnts.push_back(aPoint);

    }

    for (Integer i = 0; i < m_mesher.mesh().nbElements(); ++i)
    {

        Integer anElementId = m_mesher.mesh().elementId(i);

        CMshElement<double>& anElement = m_mesher.mesh().element(anElementId);

        if (anElement.elementType() == ET_TRIANGLE)
        {

            Integer aNodeId1 = anElement.nodeId(0);
            Integer aNodeId2 = anElement.nodeId(1);
            Integer aNodeId3 = anElement.nodeId(2);

            Integer n1 = m_mesher.mesh().nodeIndex(aNodeId1);
            Integer n2 = m_mesher.mesh().nodeIndex(aNodeId2);
            Integer n3 = m_mesher.mesh().nodeIndex(aNodeId3);

            Poly_Triangle tri(n1 + 1, n2 + 1, n3 + 1);
            tris.push_back(tri);

        }

    }

    TColgp_Array1OfPnt tab(1, pnts.size());

    for (Integer j = 0; j < pnts.size(); ++j)
        tab.SetValue(j + 1, pnts[j]);

    Poly_Array1OfTriangle trian(1, tris.size()); 

    for (Integer j = 0; j < tris.size(); ++j)
        trian.SetValue(j + 1, tris[j]);

    Handle(Poly_Triangulation) polys = new Poly_Triangulation(tab, trian);
    Handle(AIS_Triangulation) triang = new AIS_Triangulation(polys);

    m_occWindow.context()->Display(triang, Standard_True);

}

void OccManager::exportMesh(int nFormat, QString& aFileName)
{

    if (nFormat == OccTranslator::FormatGmsh)
    {

        CPdeField<double> aField;
        CPosGmsh<double> aPosGmsh;

        aField.setMesh(m_mesher.mesh());
        aField.mesh().invert();
        aField.mesh().scale(0.001);

        aPosGmsh.save(aField, aFileName.toStdString(), "");

    }
    if (nFormat == OccTranslator::FormatVtk)
    {

        CPdeField<double> aField;
        CPosVtk<double> aPosVtk;

        aField.setMesh(m_mesher.mesh());
        aPosVtk.save(aField, aFileName.toStdString());

    }
    else if (nFormat == OccTranslator::FormatSTL)
    {

        CStlUtils<double> aStlUtils(m_mesher.mesh());
        aStlUtils.save(aFileName.toStdString());

    }

}

void OccManager::selectVertices() 
{

    m_occWindow.selectVertices();

}

void OccManager::selectEdges() 
{

    m_occWindow.selectEdges();

}

void OccManager::selectWires()
{

    m_occWindow.selectWires();

}

void OccManager::selectFaces() 
{

    m_occWindow.selectFaces();

}

void OccManager::selectSolids() 
{

    m_occWindow.selectSolids();

}

void OccManager::selectAll() 
{

    m_occWindow.selectAll();

}

void OccManager::deleteSelected() 
{

    std::vector<Handle(AIS_InteractiveObject)> sObjs;

    for (m_occWindow.context()->InitCurrent(); m_occWindow.context()->MoreCurrent(); m_occWindow.context()->NextCurrent())
        sObjs.push_back(m_occWindow.context()->Current());

    for (Integer i = 0; i < sObjs.size(); ++i)
        m_occWindow.context()->Remove(sObjs[i]);

}

void OccManager::viewWireframe() 
{

    m_occWindow.viewWireframe();

}

void OccManager::viewShaded() 
{

    m_occWindow.viewShaded();

}

void OccManager::viewFront() 
{

    m_occWindow.viewFront();

}

void OccManager::viewBack() 
{

    m_occWindow.viewBack();

}

void OccManager::viewLeft() 
{

    m_occWindow.viewLeft();

}

void OccManager::viewRight() 
{

    m_occWindow.viewRight();

}

void OccManager::viewTop() 
{

    m_occWindow.viewTop();

}

void OccManager::viewBottom() 
{

    m_occWindow.viewBottom();

}

void OccManager::viewIsometric() 
{

    m_occWindow.viewIsometric();

}

void OccManager::setColor(QColor aColor)
{

    m_occWindow.setColor(aColor);

}

void OccManager::decomposeShape()
{

    try
    {

        std::vector<Handle(AIS_InteractiveObject)> sObjs;

        m_occWindow.context()->InitSelected();

        while (m_occWindow.context()->MoreSelected())
        {

            if (m_occWindow.context()->HasSelectedShape())
            {

                TopoDS_Shape aShape = m_occWindow.context()->SelectedShape();

                if (aShape.ShapeType() == TopAbs_COMPOUND || aShape.ShapeType() == TopAbs_COMPSOLID || aShape.ShapeType() == TopAbs_SOLID || aShape.ShapeType() == TopAbs_SHELL)
                {

                    for (TopExp_Explorer aSurfExp(aShape, TopAbs_FACE); aSurfExp.More(); aSurfExp.Next())
                    {

                        TopoDS_Face aFace = TopoDS::Face(aSurfExp.Current());

                        this->displayShape(aFace);

                    }

                    sObjs.push_back(m_occWindow.context()->Current());

                }

            }

            m_occWindow.context()->NextSelected();

        }

        for (Integer i = 0; i < sObjs.size(); ++i)
            m_occWindow.context()->Remove(sObjs[i]);

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

void OccManager::fuseShape()
{

    try
    {

		std::vector<Handle(AIS_InteractiveObject)> sObjs;

        m_occWindow.context()->InitSelected();

        if (m_occWindow.context()->MoreCurrent())
        {

            TopoDS_Shape aShape1 = m_occWindow.context()->SelectedShape();
			sObjs.push_back(m_occWindow.context()->Current());
			m_occWindow.context()->NextCurrent();

            if (m_occWindow.context()->MoreCurrent())
            {

                TopoDS_Shape aShape2 = m_occWindow.context()->SelectedShape();
				sObjs.push_back(m_occWindow.context()->Current());
				m_occWindow.context()->NextCurrent();

                TopoDS_Shape aJoin = BRepAlgoAPI_Fuse(aShape1, aShape2);

                if (!aJoin.IsNull())
                {
                    this->displayShape(aJoin);
                }

            }

        }

		for (Integer i = 0; i < sObjs.size(); ++i)
			m_occWindow.context()->Remove(sObjs[i]);

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

void OccManager::offsetShape(double offsetValue)
{

    try
    {

        m_occWindow.context()->InitSelected();

        while (m_occWindow.context()->MoreSelected())
        {

            if (m_occWindow.context()->HasSelectedShape())
            {
				
				TopoDS_Shape aShape = m_occWindow.context()->SelectedShape();
				
				BRepOffsetAPI_MakeOffsetShape aNewShape(aShape, offsetValue, Precision::Confusion());

				aNewShape.Build();

				this->displayShape(aNewShape.Shape());

            }

            m_occWindow.context()->NextSelected();

        }

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

void OccManager::buildShape()
{

    try
    {

        m_occWindow.context()->InitSelected();

        BRepFill_Filling aFaceBuilder;

        while (m_occWindow.context()->MoreSelected())
        {

            if (m_occWindow.context()->HasSelectedShape())
            {

                TopoDS_Shape aShape = m_occWindow.context()->SelectedShape();

                if (aShape.ShapeType() == TopAbs_EDGE)
                {

                    aFaceBuilder.Add(TopoDS::Edge(aShape), GeomAbs_C0);

                }

            }

            m_occWindow.context()->NextSelected();

        }

        aFaceBuilder.Build();

        if (aFaceBuilder.IsDone())
        {

            TopoDS_Face aFace = TopoDS::Face(aFaceBuilder.Face());

            aFace.Checked(Standard_True);

            this->displayShape(aFace);

        }

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

void OccManager::meshSurface(double meshSize)
{

    try
    {

        TopTools_DataMapOfShapeReal anEdgeMap;

        m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeOn();
        m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetInteriorStyle(Aspect_IS_EMPTY);
        m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeColor(Quantity_NOC_DARKGREEN);
        m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeLineType(Aspect_TOL_SOLID);
        m_occWindow.context()->DefaultDrawer()->ShadingAspect()->Aspect()->SetEdgeWidth(1.0);

        m_occWindow.context()->InitSelected();

        while (m_occWindow.context()->MoreSelected())
        {

            if (m_occWindow.context()->HasSelectedShape())
            {

                TopoDS_Shape aShape = m_occWindow.context()->SelectedShape();

                Integer i = 0;

                for (TopExp_Explorer aSurfExp(aShape, TopAbs_FACE); aSurfExp.More(); aSurfExp.Next())
                {

                    TopoDS_Face aFace = TopoDS::Face(aSurfExp.Current());

                    std::cout << "aFace = " << i << std::endl;
                    i++;

                    std::vector<gp_Pnt> pnts;
                    std::vector<Poly_Triangle> tris;

                    bool res = m_mesher.meshSurface(aFace, pnts, tris, meshSize, 1E-2);

                    if (res)
                    {
                        TColgp_Array1OfPnt tab(1, pnts.size());

                        for (Integer j = 0; j < pnts.size(); ++j)
                            tab.SetValue(j + 1, pnts[j]);

                        Poly_Array1OfTriangle trian(1, tris.size()); 

                        for (Integer j = 0; j < tris.size(); ++j)
                            trian.SetValue(j + 1, tris[j]);

                        Handle(Poly_Triangulation) polys = new Poly_Triangulation(tab, trian);
                        Handle(AIS_Triangulation) triang = new AIS_Triangulation(polys);

                        m_occWindow.context()->Display(triang, Standard_True);

                    }

                }

            }

            m_occWindow.context()->NextSelected();

        }

        m_mesher.mesh().mergeNodes(1E-2);
        m_mesher.mesh().renumber();

        Integer n = 0;

        for (Integer i = 0; i < m_mesher.mesh().nbElements(); ++i)
        {

            Integer anElementId = m_mesher.mesh().elementId(i);

            if (m_mesher.mesh().element(anElementId).elementType() == ET_TRIANGLE ||
                m_mesher.mesh().element(anElementId).elementType() == ET_QUADRILATERAL)
                n++;

        }

        std::cout << "Number of surface elements: " << n << std::endl;

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

void OccManager::meshVolume(double meshSize)
{

    try
    {

        m_mesher.meshVolume(meshSize, 1E-3);

        Integer n = 0;

        for (Integer i = 0; i < m_mesher.mesh().nbElements(); ++i)
        {

            Integer anElementId = m_mesher.mesh().elementId(i);

            if (m_mesher.mesh().element(anElementId).elementType() == ET_TETRAHEDRON)
                n++;
        }

        std::cout << "Number of volume elements: " << n << std::endl;

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

void OccManager::getProperties(QList<QString>& sPropertyNames, QList<QVariant>& sPropertyValues)
{

    try
    {

        m_occWindow.context()->InitSelected();

        while (m_occWindow.context()->MoreSelected())
        {

            if (m_occWindow.context()->HasSelectedShape())
            {

                TopoDS_Shape aShape = m_occWindow.context()->SelectedShape();

                GProp_GProps sProps;

                if (aShape.ShapeType() == TopAbs_EDGE || aShape.ShapeType() == TopAbs_WIRE)
                {
                    if (aShape.ShapeType() == TopAbs_EDGE)
                    {
                        sPropertyNames.push_back(tr("Edge"));
                        sPropertyValues.push_back("");
                    }

                    BRepGProp::LinearProperties(aShape, sProps);
                    sPropertyNames.push_back(tr("Length"));
                    sPropertyValues.push_back(sProps.Mass());
                }
                else if (aShape.ShapeType() == TopAbs_FACE || aShape.ShapeType() == TopAbs_SHELL)
                {
                    BRepGProp::SurfaceProperties(aShape, sProps);
                    sPropertyNames.push_back(tr("Area"));
                    sPropertyValues.push_back(sProps.Mass());
                }

            }

            break;

        }

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << tr("Unknown exception:") << QString(this->metaObject()->className());
    }

}

void OccManager::meshStitch()
{

    m_mesher.stitchMesh();

}

void OccManager::selectShapes(QModelIndexList selected)
{

    m_occWindow.context()->ClearSelected();

    m_occWindow.context()->InitSelected();

    for(QList<QModelIndex>::const_iterator iter = selected.constBegin(); iter != selected.constEnd(); ++iter)
    {

        TopoDS_Shape aShape = m_shapeMap[(*iter)];

        m_occWindow.context()->AddOrRemoveSelected(aShape);

    }

}

void OccManager::treeViewSelectionChanged(const QItemSelection& selected, const QItemSelection& deselected)
{

    try
    {

        this->selectShapes(selected.indexes());

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
