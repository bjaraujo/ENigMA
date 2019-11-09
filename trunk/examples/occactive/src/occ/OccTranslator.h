#pragma once

#include <iostream>
#include <string>

#include <TopTools_HSequenceOfShape.hxx>

class OccTranslator
{
public:

    enum {FormatBREP, FormatIGES, FormatSTEP, FormatCSFDB, FormatVRML, FormatSTL, FormatGmsh, FormatVtk};

    OccTranslator();
    ~OccTranslator();

    std::string info() const;

    Handle(TopTools_HSequenceOfShape) importModel(const int nFormat, const std::string file);
    bool exportModel(const int format, const std::string& file, const Handle(TopTools_HSequenceOfShape)& shapes);

private:

    Handle(TopTools_HSequenceOfShape) importBREP(const std::string& fileName);
    Handle(TopTools_HSequenceOfShape) importIGES(const std::string& fileName);
    Handle(TopTools_HSequenceOfShape) importSTEP(const std::string& fileName);
    Handle(TopTools_HSequenceOfShape) importCSFDB(const std::string& fileName);

    bool exportBREP(const std::string&, const Handle(TopTools_HSequenceOfShape)&);
    bool exportIGES(const std::string&, const Handle(TopTools_HSequenceOfShape)&);
    bool exportSTEP(const std::string&, const Handle(TopTools_HSequenceOfShape)&);
    bool exportCSFDB(const std::string&, const Handle(TopTools_HSequenceOfShape)&);
    bool exportSTL(const std::string&, const Handle(TopTools_HSequenceOfShape)&);
    bool exportVRML(const std::string&, const Handle(TopTools_HSequenceOfShape)&);

    bool checkFacetedBrep(const Handle(TopTools_HSequenceOfShape)&);

protected:
    std::string        m_info;

};


