#pragma once

#include <iostream>
#include <string>

#include <TopTools_DataMapOfShapeReal.hxx>

#include "MshCheckMesh.hpp"
#include "MshTriangleMesher.hpp"
#include "MshQuadrilateralMesher.hpp"
#include "MshTetrahedronMesher.hpp"
#include "StlUtils.hpp"
#include "PdeEquation.hpp"
#include "PosGmsh.hpp"
#include "PosVtk.hpp"
#include "PosQuickMesh.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::pde;
using namespace ENigMA::post;
using namespace ENigMA::stl;

class TopoDS_Face;
class TopoDS_Edge;
class gp_Pnt;
class gp_Pnt2d;
class Poly_Triangle;
class TopTools_IndexedMapOfShape;
	 
class OccMesher
{
public:
    OccMesher();
    ~OccMesher();

    void reset();

    bool meshSurface(const TopoDS_Face& aFace, std::vector<gp_Pnt>& pnts, std::vector<Poly_Triangle>& tris, const double h, const double tolerance);
    bool meshVolume(double h, double tolerance);

    void stitchMesh();

    CMshMesh<double>& mesh() { return m_mesh; }

private:
    CMshMesh<double> m_mesh;

    TopTools_DataMapOfShapeReal m_edgeMap;
    TopTools_DataMapOfShapeReal m_faceMap;

private:

	void discretizeEdge(CMshMesh<double>& anEdgeMesh, const TopoDS_Face& aFace, const TopoDS_Edge& anEdge, CGeoRtree<double>& aRtree, Integer& aLastFaceId, const double aScaleX, const double aScaleY, const double offsetZ, const bool bReverse, const double h, const double tolerance);

    bool getSurfaceUVFromPoint(const TopoDS_Face& aFace, gp_Pnt aPoint, gp_Pnt2d& a2DPoint, const double tolerance);
    bool getPointFromSurfaceUV(const TopoDS_Face& aFace, gp_Pnt2d a2DPoint, gp_Pnt& aPoint);

	Integer addPointToMesh(CMshMesh<double>& anEdgeMesh, gp_Pnt& aPoint, CGeoRtree<double>& aRtree, const double tolerance);

    void computeScaleOnFace(const TopoDS_Face& aFace, double& scalex, double& scaley);
    
};


