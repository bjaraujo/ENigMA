#if !defined(CSFDB)
#error CSFDB precompiler directive is mandatory for CasCade 
#endif

#include <AIS_Shape.hxx>
#include <AIS_InteractiveObject.hxx>

#include <FSD_File.hxx>

#include <ShapeSchema.hxx>
#include <Storage_Data.hxx>
#include <Storage_Root.hxx>
#include <Storage_HSeqOfRoot.hxx>
#include <PTopoDS_HShape.hxx>
#include <PTColStd_PersistentTransientMap.hxx>
#include <PTColStd_TransientPersistentMap.hxx>

#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <IGESControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <Interface_Static.hxx>

#include <StlAPI_Writer.hxx>
#include <VrmlAPI_Writer.hxx>

#include <MgtBRep.hxx>
#include <MgtBRep_TriangleMode.hxx>

#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_HSequenceOfShape.hxx>

#include <Geom_Line.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>

#include <Standard_ErrorHandler.hxx>
#include <Standard_CString.hxx>

#include "OccTranslator.h"

OccTranslator::OccTranslator()
{

}

OccTranslator::~OccTranslator()
{

}

Handle(TopTools_HSequenceOfShape) OccTranslator::importModel(const int nFormat, const std::string file)
{

    Handle(TopTools_HSequenceOfShape) shapes;

    try 
    {
        switch (nFormat)
        {
        case FormatBREP:
            shapes = importBREP(file);
            break;
        case FormatIGES:
            shapes = importIGES(file);
            break;
        case FormatSTEP:
            shapes = importSTEP(file);
            break;
        case FormatCSFDB:
            shapes = importCSFDB(file);
            break;
        }
    } 
    catch (Standard_Failure) 
    {
        shapes.Nullify();
    }

    return shapes;

}

bool OccTranslator::exportModel(const int format, const std::string& file, const Handle(TopTools_HSequenceOfShape)& shapes)
{

    bool status;
    try 
    {
        switch (format)
        {
        case FormatBREP:
            status = exportBREP(file, shapes);
            break;
        case FormatIGES:
            status = exportIGES(file, shapes);
            break;
        case FormatSTEP:
            status = exportSTEP(file, shapes);
            break;
        case FormatCSFDB:
            status = exportCSFDB(file, shapes);
            break;
        case FormatSTL:
            status = exportSTL(file, shapes);
            break;
        case FormatVRML:
            status = exportVRML(file, shapes);
            break;
        }
    } 
    catch (Standard_Failure) 
    {
        status = false;
    }

    return status;
}

Handle(TopTools_HSequenceOfShape) OccTranslator::importBREP(const std::string& file)
{
    Handle(TopTools_HSequenceOfShape) aSequence;
    TopoDS_Shape aShape;
    BRep_Builder aBuilder;

    Standard_Boolean result = BRepTools::Read(aShape, (Standard_CString)file.c_str(), aBuilder);
    if (result)
    {
        aSequence = new TopTools_HSequenceOfShape();
        aSequence->Append(aShape);
    }
    return aSequence;
}

Handle(TopTools_HSequenceOfShape) OccTranslator::importIGES(const std::string& file)
{
    Handle(TopTools_HSequenceOfShape) aSequence;
    IGESControl_Reader Reader;
    int status = Reader.ReadFile((Standard_CString)file.c_str());

    if (status == IFSelect_RetDone)
    {
        aSequence = new TopTools_HSequenceOfShape();
        Reader.TransferRoots();
        TopoDS_Shape aShape = Reader.OneShape();
        aSequence->Append(aShape);
    }
    return aSequence;
}

Handle(TopTools_HSequenceOfShape) OccTranslator::importSTEP(const std::string& file)
{
    Handle(TopTools_HSequenceOfShape) aSequence;

    STEPControl_Reader aReader;
    IFSelect_ReturnStatus status = aReader.ReadFile((Standard_CString)file.c_str());

    if (status == IFSelect_RetDone)
    {
        //Interface_TraceFile::SetDefault();
        //bool failsonly = false;
        //aReader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity);

        int nbr = aReader.NbRootsForTransfer();
        //aReader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity);
        for (Standard_Integer n = 1; n <= nbr; n++)
        {
            aReader.TransferRoot(n);
            int nbs = aReader.NbShapes();
            if (nbs > 0)
            {
                aSequence = new TopTools_HSequenceOfShape();
                for (int i = 1; i <= nbs; ++i)
                {
                    TopoDS_Shape aShape = aReader.Shape(i);
                    aSequence->Append(aShape);
                }
            }
        }
    }
    return aSequence;
}

Handle(TopTools_HSequenceOfShape) OccTranslator::importCSFDB(const std::string& file)
{
    Handle(TopTools_HSequenceOfShape) aSequence;

    // Check file type
    if (FSD_File::IsGoodFileType((Standard_CString)file.c_str()) != Storage_VSOk)
        return aSequence;

    static FSD_File fileDriver;
    TCollection_AsciiString aName((Standard_CString)file.c_str());
    if (fileDriver.Open(aName, Storage_VSRead) != Storage_VSOk)
        return aSequence;

    Handle(ShapeSchema) schema = new ShapeSchema();
    Handle(Storage_Data) data  = schema->Read(fileDriver);
    if (data->ErrorStatus() != Storage_VSOk)
        return aSequence;

    fileDriver.Close();

    aSequence = new TopTools_HSequenceOfShape();
    Handle(Storage_HSeqOfRoot) roots = data->Roots();
    for (int i = 1; i <= roots->Length() ; ++i)
    {
        Handle(Storage_Root) r = roots->Value(i);
        Handle(Standard_Persistent) p = r->Object();
        Handle(PTopoDS_HShape) aPShape = Handle(PTopoDS_HShape)::DownCast(p);
        if (!aPShape.IsNull())
        {
            PTColStd_PersistentTransientMap aMap;
            TopoDS_Shape aTShape;
            MgtBRep::Translate(aPShape, aMap, aTShape, MgtBRep_WithTriangle);
            aSequence->Append(aTShape);
        }
    }

    return aSequence;
}

// ----------------------------- Export functionality -----------------------------

bool OccTranslator::exportBREP(const std::string& file, const Handle(TopTools_HSequenceOfShape)& shapes)
{

    if (shapes.IsNull() || shapes->IsEmpty())
        return false;

    TopoDS_Shape aShape = shapes->Value(1);
    return BRepTools::Write(aShape, (Standard_CString)file.c_str()) == Standard_True ? true : false;

}

bool OccTranslator::exportIGES(const std::string& file, const Handle(TopTools_HSequenceOfShape)& shapes)
{
    if (shapes.IsNull() || shapes->IsEmpty())
        return false;

    IGESControl_Controller::Init();
    IGESControl_Writer writer(Interface_Static::CVal("XSTEP.iges.unit"),
        Interface_Static::IVal("XSTEP.iges.writebrep.mode"));

    for (int i = 1; i <= shapes->Length(); ++i)
        writer.AddShape (shapes->Value(i));
    writer.ComputeModel();

    return writer.Write((Standard_CString)file.c_str()) == Standard_True ? true : false;

}

bool OccTranslator::exportSTEP(const std::string& file, const Handle(TopTools_HSequenceOfShape)& shapes)
{

    if (shapes.IsNull() || shapes->IsEmpty())
        return false;

    IFSelect_ReturnStatus status;

    STEPControl_Writer writer;
    for (int i = 1; i <= shapes->Length(); ++i)
    {
        status = writer.Transfer(shapes->Value(i), STEPControl_AsIs);
        if (status != IFSelect_RetDone)
            return false;
    }

    status = writer.Write((Standard_CString)file.c_str());

    switch (status)
    {
    case IFSelect_RetError:
        m_info = "Incorrect data.";
        break;
    case IFSelect_RetFail:
        m_info = "Write fail.";
        break;
    case IFSelect_RetVoid:
        m_info = "Nothing to transfer.";
        break;
    }
    return status == IFSelect_RetDone;

}

bool OccTranslator::exportCSFDB(const std::string& file, const Handle(TopTools_HSequenceOfShape)& shapes)
{

    if (shapes.IsNull() || shapes->IsEmpty())
        return false;

    static FSD_File fileDriver;

    Handle(ShapeSchema) schema = new ShapeSchema();
    Handle(Storage_Data) data  = new Storage_Data();
    data->ClearErrorStatus();

    /*
    data->SetApplicationName(TCollection_ExtendedString("Sample Import / Export"));
    data->SetApplicationVersion("1");
    data->SetDataType(TCollection_ExtendedString("Shapes"));
    data->AddToUserInfo("Storing a persistent set of shapes in a flat file");
    data->AddToComments(TCollection_ExtendedString("Application is based on CasCade 5.0 Professional"));
    */

    if (fileDriver.Open((Standard_CString)file.c_str(), Storage_VSWrite) != Storage_VSOk)
    {
        m_info = "Can't save file: " + file;
        return false;
    }

    PTColStd_TransientPersistentMap aMap;
    for (int i = 1; i <= shapes->Length(); ++i)
    {
        TopoDS_Shape aShape = shapes->Value(i);
        if (aShape.IsNull())
        {
            m_info = "Some shapes are invalid.";
            return false;
        }

        Handle(PTopoDS_HShape) pshape = MgtBRep::Translate(aShape, aMap, MgtBRep_WithTriangle);
        TCollection_AsciiString objName = TCollection_AsciiString("Object_") + TCollection_AsciiString(i);
        data->AddRoot(objName, pshape);
    }

    schema->Write(fileDriver, data);
    fileDriver.Close();

    if (data->ErrorStatus() != Storage_VSOk)
    {
        m_info = "Can't store persistent data.";
        return false;
    } 

    return true;
}

bool OccTranslator::exportSTL(const std::string& file, const Handle(TopTools_HSequenceOfShape)& shapes)
{
    if (shapes.IsNull() || shapes->IsEmpty())
        return false;

    TopoDS_Compound res;
    BRep_Builder builder;
    builder.MakeCompound(res);

    for (int i = 1; i <= shapes->Length(); ++i)
    {
        TopoDS_Shape aShape = shapes->Value(i);
        if (aShape.IsNull()) 
        {
            m_info = "Some shapes are invalid.";
            return false;
        }
        builder.Add(res, aShape);
    }

    StlAPI_Writer writer;
    writer.Write(res, (Standard_CString)file.c_str());

    return true;
}

bool OccTranslator::exportVRML(const std::string& file, const Handle(TopTools_HSequenceOfShape)& shapes)
{
    if (shapes.IsNull() || shapes->IsEmpty())
        return false;

    TopoDS_Compound res;
    BRep_Builder builder;
    builder.MakeCompound(res);

    for (int i = 1; i <= shapes->Length(); ++i)
    {
        TopoDS_Shape aShape = shapes->Value(i);
        if (aShape.IsNull())
        {
            m_info = "Some shapes are invalid.";
            return false;
        }
        builder.Add(res, aShape);
    }

    VrmlAPI_Writer writer;
    writer.Write(res, (Standard_CString)file.c_str());

    return true;
}

bool OccTranslator::checkFacetedBrep(const Handle(TopTools_HSequenceOfShape)& shapes)
{
    bool err = false;
    for (int i = 1; i <= shapes->Length(); ++i)
    {
        TopoDS_Shape aShape = shapes->Value(i);
        for (TopExp_Explorer fexp(aShape, TopAbs_FACE); fexp.More() && !err; fexp.Next())
        {
            Handle(Geom_Surface) surface = BRep_Tool::Surface(TopoDS::Face(fexp.Current()));
            if (!surface->IsKind(STANDARD_TYPE(Geom_Plane)))
                err = true;
        }
        for (TopExp_Explorer eexp(aShape, TopAbs_EDGE); eexp.More() && !err; eexp.Next())
        {
            Standard_Real fd, ld;
            Handle(Geom_Curve) curve = BRep_Tool::Curve(TopoDS::Edge(eexp.Current()), fd, ld);
            if (!curve->IsKind(STANDARD_TYPE(Geom_Line)))
                err = true;
        }
    }
    return !err;
}



