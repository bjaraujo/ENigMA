
#include "OccMesher.h"

#include <Precision.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Curve.hxx>
#include <Geom2d_Curve.hxx>

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ShapeAnalysis_Curve.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <ShapeExtend_WireData.hxx>
#include <BRepTools.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepTools_WireExplorer.hxx>

#include <Geom_CartesianPoint.hxx>
#include <GeomLProp_CLProps.hxx>
#include <GeomLProp_SLProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>

#include <TShort_HArray1OfShortReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <Poly_Array1OfTriangle.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_Polygon3D.hxx>
#include <Poly_PolygonOnTriangulation.hxx>

#include <AIS_Point.hxx>
#include <AIS_Line.hxx>
#include <AIS_Triangulation.hxx>

#include <Prs3d_Drawer.hxx>
#include <Prs3d_ShadingAspect.hxx>
#include <Graphic3d_AspectFillArea3d.hxx>

#include <QElapsedTimer>

#ifdef USE_NETGEN
namespace nglib {
#include "nglib.h"
}

using namespace nglib;
#endif

OccMesher::OccMesher()
{

}

OccMesher::~OccMesher()
{


}

void OccMesher::reset()
{

    m_mesh.reset();
    m_edgeMap.Clear();
    m_faceMap.Clear();

}

void OccMesher::computeScaleOnFace(const TopoDS_Face& aFace, double& aScaleX, double& aScaleY)
{

    TopoDS_Wire W = BRepTools::OuterWire(aFace);

    double xmin = +std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = +std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    int nbp = 23;
    aScaleX = 1.0;
    aScaleY = 1.0;

    TopExp_Explorer wexp(W, TopAbs_EDGE);

    for ( ; wexp.More(); wexp.Next())
    {
        const TopoDS_Edge & E = TopoDS::Edge( wexp.Current() );
        double f, l;
        Handle(Geom2d_Curve) C2d = BRep_Tool::CurveOnSurface(E, aFace, f, l);
        if ( C2d.IsNull() ) continue;
        double du = (l - f) / double (nbp);
        for (int i = 0; i <= nbp; ++i)
        {
            double param = f + double (i) * du;
            gp_Pnt2d p = C2d->Value(param);
            if (p.X() < xmin)
                xmin = p.X();
            if (p.X() > xmax)
                xmax = p.X();
            if (p.Y() < ymin)
                ymin = p.Y();
            if (p.Y() > ymax)
                ymax = p.Y();
        }

    }

    double xmoy = (xmax + xmin) / 2.;
    double ymoy = (ymax + ymin) / 2.;
    double xsize = xmax - xmin;
    double ysize = ymax - ymin;

    TopLoc_Location L;
    Handle(Geom_Surface) S = BRep_Tool::Surface(aFace, L);       // 3D surface

    double length_x = 0.0;
    double length_y = 0.0;

    gp_Pnt PX0 = S->Value(xmin, ymoy);
    gp_Pnt PY0 = S->Value(xmoy, ymin);

    double dx = xsize / double (nbp);
    double dy = ysize / double (nbp);

    for (int i = 1; i <= nbp; ++i)
    {
        double x = xmin + double (i) * dx;
        gp_Pnt PX = S->Value(x, ymoy);
        double y = ymin + double (i) * dy;
        gp_Pnt PY = S->Value(xmoy, y);
        length_x += PX.Distance(PX0);
        length_y += PY.Distance(PY0);
        PX0 = PX;
        PY0 = PY;
    }

    aScaleX = length_x / xsize;
    aScaleY = length_y / ysize;

    double xyratio = xsize*aScaleX/(ysize*aScaleY);
    const double maxratio = 1.e2;

    if (xyratio > maxratio) 
    {
        aScaleY *= xyratio / maxratio;
    }
    else if (xyratio < 1./maxratio) 
    {
        aScaleX *= 1 / xyratio / maxratio;
    }

}

bool OccMesher::getSurfaceUVFromPoint(const TopoDS_Face& aFace, gp_Pnt aPoint, gp_Pnt2d& a2DPoint, const double tolerance)
{

    Handle(Geom_Surface) aSurf = BRep_Tool::Surface(aFace);
    Handle(ShapeAnalysis_Surface) aSurfAna = new ShapeAnalysis_Surface(aSurf);
    a2DPoint = aSurfAna->ValueOfUV(aPoint, tolerance);

    double u, v;
    a2DPoint.Coord(u, v);

    double umin, umax, vmin, vmax;

    BRepTools::UVBounds(aFace, umin, umax, vmin, vmax); 

    if (u < umin || u > umax || v < vmin || v > vmax)
    {
        cout << "u =" << u << " v = " << v << endl;
        cout << "Bounds: u = [" << umin << ", " << umax << "] v = [" << vmin << ", " << vmax << "]" << endl;
        cout << "Error: out of bounds! " << endl;
        return false;
    }
    else
        return true;

}

bool OccMesher::getPointFromSurfaceUV(const TopoDS_Face& aFace, gp_Pnt2d a2DPoint, gp_Pnt& aPoint)
{

    Handle_Geom_Surface aSurf = BRep_Tool::Surface(aFace);

    double u, v;
    a2DPoint.Coord(u, v);

    aPoint = aSurf->Value(u, v);

    return true;

}

Integer OccMesher::addPointToMesh(CMshMesh<double>& aMesh, gp_Pnt& aPoint, CGeoRtree<double>& aRtree, double tolerance)
{

    CMshNode<double> aNode(aPoint.X(), aPoint.Y(), aPoint.Z());

    CGeoBoundingBox<double> aBoundingBox;

    aBoundingBox.addCoordinate(aNode);
    aBoundingBox.grow(tolerance);

    std::vector<Integer> sNodes;

    aRtree.find(sNodes, aBoundingBox, tolerance);

    if (sNodes.size() == 0)
    {

        // Create node
        Integer aNodeId = aMesh.nextNodeId();
        aMesh.addNode(aNodeId, aNode);
        aRtree.addGeometricObject(aNodeId, aBoundingBox);

        return aNodeId;

    }
    else
    {

        for (Integer i = 0; i < sNodes.size(); ++i)
        {

            double d = (aNode - aMesh.node(sNodes[i])).norm();

            if (d < tolerance)
                return sNodes[i];

        }

        return sNodes[0];

    }

}

void OccMesher::discretizeEdge(CMshMesh<double>& anEdgeMesh, const TopoDS_Face& aFace, const TopoDS_Edge& anEdge, CGeoRtree<double>& aRtree, Integer& aLastFaceId, const double aScaleX, const double aScaleY, const double offsetZ, const bool bReverse, const double h, const double tolerance)
{

    bool bReverseEdge = (aFace.Orientation() == TopAbs_REVERSED) ^ (anEdge.Orientation() == TopAbs_REVERSED); 

    if (bReverse)
        bReverseEdge = !bReverseEdge;

    double s0, s1;

    double aCurvature;

    Handle(Geom_Curve) aCurve = BRep_Tool::Curve(anEdge, s0, s1);

    if (!aCurve.IsNull())
    {

        GeomLProp_CLProps props(aCurve, (s0 + s1)*0.5, 2, Precision::Confusion());
        aCurvature = props.Curvature();

    }
    else
        aCurvature = 0.0;

    GProp_GProps System;
    BRepGProp::LinearProperties(anEdge, System);
    Standard_Real Length = System.Mass();

    Handle(Geom2d_Curve) a2DCurve = BRep_Tool::CurveOnSurface(anEdge, aFace, s0, s1);

    gp_Pnt2d a2DPoint1, a2DPoint2;

    if (bReverseEdge)
        a2DPoint1 = a2DCurve->Value(s1);
    else
        a2DPoint1 = a2DCurve->Value(s0);

    a2DPoint1.SetX(a2DPoint1.X() * aScaleX);
    a2DPoint1.SetY(a2DPoint1.Y() * aScaleY);

	Integer aFirstNodeId = this->addPointToMesh(anEdgeMesh, gp_Pnt(a2DPoint1.X(), a2DPoint1.Y(), offsetZ), aRtree, tolerance);

    int nDiv = static_cast<int> (ceil(Length / h));

    const double pi = std::acos(-1.0);

    BRepAdaptor_Curve anAdaptorCurve(anEdge);

    if (anAdaptorCurve.GetType() == GeomAbs_Circle)
    {

        //cout << "Circular edge " << endl;
        nDiv = std::max((int) (fabs(s1 - s0) / pi * 3), nDiv);

    }
    else
    {

        if (aCurvature > 0.0 && nDiv < 2) 
            nDiv = 2;

    }

    double d = (s1 - s0) / nDiv;

    Integer aPrevNodeId = aFirstNodeId;

    for (int i = 0; i < nDiv - 1; ++i)
    {

        gp_Pnt2d a2DPoint;

        if (bReverseEdge)
            a2DPoint = a2DCurve->Value(s1 - (i + 1) * d);
        else
            a2DPoint = a2DCurve->Value(s0 + (i + 1) * d);

        a2DPoint.SetX(a2DPoint.X() * aScaleX);
        a2DPoint.SetY(a2DPoint.Y() * aScaleY);

		Integer aNodeId = this->addPointToMesh(anEdgeMesh, gp_Pnt(a2DPoint.X(), a2DPoint.Y(), offsetZ), aRtree, tolerance);

        if (aPrevNodeId != aNodeId)
        {

            CMshElement<double> anElement(ET_BEAM);

            anElement.addNodeId(aPrevNodeId);
            anElement.addNodeId(aNodeId);

            Integer anElementId = anEdgeMesh.nextElementId();

            std::vector<CMshFace<double> > sFaces;
            anElement.generateFaces(sFaces);

            for (Integer j = 0; j < sFaces.size(); ++j)
            {

                CMshFace<double> aFace = sFaces[j];

                aFace.setElementId(anElementId);

                Integer aFaceId = anEdgeMesh.nextFaceId();

                if (j == 0)
                {
                    if (aFaceId > 0)
                        aFace.setPairFaceId(aFaceId - 1);
                }
                else
                    aFace.setPairFaceId(aFaceId + 1);

                anEdgeMesh.addFace(aFaceId, aFace);
                anElement.addFaceId(aFaceId);

            }

            anEdgeMesh.addElement(anElementId, anElement);

        }

        aPrevNodeId = aNodeId;

    }

    if (bReverseEdge)
        a2DPoint2 = a2DCurve->Value(s0);
    else
        a2DPoint2 = a2DCurve->Value(s1);

    a2DPoint2.SetX(a2DPoint2.X() * aScaleX);
    a2DPoint2.SetY(a2DPoint2.Y() * aScaleY);

	Integer aLastNodeId = this->addPointToMesh(anEdgeMesh, gp_Pnt(a2DPoint2.X(), a2DPoint2.Y(), offsetZ), aRtree, tolerance);

    CMshElement<double> anElement(ET_BEAM);

    anElement.addNodeId(aPrevNodeId);
    anElement.addNodeId(aLastNodeId);

    Integer anElementId = anEdgeMesh.nextElementId();

    std::vector<CMshFace<double> > sFaces;
    anElement.generateFaces(sFaces);

    for (Integer j = 0; j < sFaces.size(); ++j)
    {

        CMshFace<double> aFace = sFaces[j];

        aFace.setElementId(anElementId);

        Integer aFaceId = anEdgeMesh.nextFaceId();

        if (j == 0)
        {
            if (aFaceId > 0)
                aFace.setPairFaceId(aFaceId - 1);
        }
        else
            aFace.setPairFaceId(aFaceId + 1);

        anEdgeMesh.addFace(aFaceId, aFace);
        anElement.addFaceId(aFaceId);

        aLastFaceId = aFaceId;

    }

    anEdgeMesh.addElement(anElementId, anElement);

}

bool OccMesher::meshSurface(const TopoDS_Face& aFace, std::vector<gp_Pnt>& pnts, std::vector<Poly_Triangle>& tris, double h, double tolerance)
{

    if (m_faceMap.IsBound(aFace))
    {
        std::cout << "Face already meshed. Skipping..." << std::endl;
        return false;
    }

    CMshMesh<double> anEdgeMesh;

    CGeoRtree<double> aRtree;

    double aScaleX, aScaleY;

    this->computeScaleOnFace(aFace, aScaleX, aScaleY);

	TopTools_IndexedMapOfShape sEdges;
	TopExp::MapShapes(aFace, TopAbs_EDGE, sEdges);

	// Internal boundaries
	for (TopExp_Explorer aWireExp(aFace, TopAbs_WIRE); aWireExp.More(); aWireExp.Next())
	{

		TopoDS_Wire aWire = TopoDS::Wire(aWireExp.Current());

		Integer aFirstFaceId = anEdgeMesh.nextFaceId();
		Integer aLastFaceId = aFirstFaceId;

		for (TopExp_Explorer anEdgeExp(aWire, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next())
		{

			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExp.Current());

			if (anEdge.Orientation() == TopAbs_INTERNAL)
			{
				if (m_edgeMap.IsBound(anEdge))
					h = m_edgeMap(anEdge);
				else
					m_edgeMap.Bind(anEdge, h);

				this->discretizeEdge(anEdgeMesh, aFace, anEdge, aRtree, aLastFaceId, aScaleX, aScaleY, 0.0, false, h, tolerance);
			}

		}

		for (TopExp_Explorer anEdgeExp(aWire, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next())
		{

			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExp.Current());

			if (anEdge.Orientation() == TopAbs_INTERNAL)
			{
				if (m_edgeMap.IsBound(anEdge))
					h = m_edgeMap(anEdge);
				else
					m_edgeMap.Bind(anEdge, h);

				this->discretizeEdge(anEdgeMesh, aFace, anEdge, aRtree, aLastFaceId, aScaleX, aScaleY, 0.0, true, h, tolerance);
			}

		}

		// Close
		anEdgeMesh.face(aFirstFaceId).setPairFaceId(aLastFaceId);
		anEdgeMesh.face(aLastFaceId).setPairFaceId(aFirstFaceId);

	}

	// Outer boundaries
    for (TopExp_Explorer aWireExp(aFace, TopAbs_WIRE); aWireExp.More(); aWireExp.Next())
    {

        TopoDS_Wire aWire = TopoDS::Wire(aWireExp.Current());

        Integer aFirstFaceId = anEdgeMesh.nextFaceId();
        Integer aLastFaceId = aFirstFaceId;

		for (TopExp_Explorer anEdgeExp(aWire, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next())
		{

			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExp.Current());

			if (anEdge.Orientation() != TopAbs_INTERNAL)
			{

				if (m_edgeMap.IsBound(anEdge))
					h = m_edgeMap(anEdge);
				else
					m_edgeMap.Bind(anEdge, h);

				this->discretizeEdge(anEdgeMesh, aFace, anEdge, aRtree, aLastFaceId, aScaleX, aScaleY, 0.0, false, h, tolerance);
			}
			
        }

        // Close
        anEdgeMesh.face(aFirstFaceId).setPairFaceId(aLastFaceId);
        anEdgeMesh.face(aLastFaceId).setPairFaceId(aFirstFaceId);

    }

	m_faceMap.Bind(aFace, h);

    QElapsedTimer timer;

    CMshTriangleMesher<double> aMesher;

    timer.restart();

	std::vector<CGeoCoordinate<double> > sInteriorPoints;

	aMesher.generate(anEdgeMesh, 99999, sInteriorPoints, h, 0.1, tolerance*tolerance);
 
	std::cout << "Elapsed time: " << timer.elapsed() / 1000.0 << " seconds" << std::endl;

	aMesher.applyFixedBoundary(anEdgeMesh, tolerance*tolerance);

    for (Integer i = 0; i < 2; ++i)
    {
		aMesher.flipEdges(tolerance*tolerance);
		aMesher.applyFixedBoundary(anEdgeMesh, tolerance*tolerance);
		aMesher.relaxNodes(tolerance*tolerance);
    }

	CMshMesh<double> aSurfaceMesh = aMesher.mesh();

    for (Integer i = 0; i < aSurfaceMesh.nbNodes(); ++i)
    {

        Integer aNodeId = aSurfaceMesh.nodeId(i);

        CMshNode<double>& aNode = aSurfaceMesh.node(aNodeId);

        gp_Pnt2d a2DPoint(aNode.x() / aScaleX, aNode.y() / aScaleY);

        gp_Pnt aPoint;

        this->getPointFromSurfaceUV(aFace, a2DPoint, aPoint);

        pnts.push_back(aPoint);

        aNode.x() = aPoint.X();
        aNode.y() = aPoint.Y();
        aNode.z() = aPoint.Z();

    }

    for (Integer i = 0; i < aSurfaceMesh.nbElements(); ++i)
    {

        Integer anElementId = aSurfaceMesh.elementId(i);

        CMshElement<double>& anElement = aSurfaceMesh.element(anElementId);

        if (aFace.Orientation() == TopAbs_REVERSED)
            anElement.invert();

        if (anElement.elementType() == ET_TRIANGLE)
        {

            Integer aNodeId1 = anElement.nodeId(0);
            Integer aNodeId2 = anElement.nodeId(1);
            Integer aNodeId3 = anElement.nodeId(2);

            Integer n1 = aSurfaceMesh.nodeIndex(aNodeId1);
            Integer n2 = aSurfaceMesh.nodeIndex(aNodeId2);
            Integer n3 = aSurfaceMesh.nodeIndex(aNodeId3);

            Poly_Triangle tri(n1 + 1, n2 + 1, n3 + 1);
            tris.push_back(tri);

        }

        if (anElement.elementType() == ET_QUADRILATERAL)
        {

            Integer aNodeId1 = anElement.nodeId(0);
            Integer aNodeId2 = anElement.nodeId(1);
            Integer aNodeId3 = anElement.nodeId(2);
            Integer aNodeId4 = anElement.nodeId(3);

            Integer n1 = aSurfaceMesh.nodeIndex(aNodeId1);
            Integer n2 = aSurfaceMesh.nodeIndex(aNodeId2);
            Integer n3 = aSurfaceMesh.nodeIndex(aNodeId3);
            Integer n4 = aSurfaceMesh.nodeIndex(aNodeId4);

            Poly_Triangle tri1(n1 + 1, n2 + 1, n3 + 1);
            Poly_Triangle tri2(n1 + 1, n3 + 1, n4 + 1);
            tris.push_back(tri1);
            tris.push_back(tri2);

        }

    }

    m_mesh.addMesh(aSurfaceMesh);

    m_mesh.renumber();

    double minq, maxq, aveq;

    m_mesh.meshQuality(minq, maxq, aveq);
    std::cout << "Mesh quality: " << aveq << " [" << minq << ", " << maxq << "]" << std::endl;

    return true;

}

bool OccMesher::meshVolume(double h, double tolerance)
{

    m_mesh.mergeNodes(tolerance);

    #ifdef USE_NETGEN
    Ng_Init();

    double pnt[3];
    int tri[3];
    int tet[4];

    Ng_Mesh* mesh = Ng_NewMesh();

    for (Integer i = 0; i < m_mesh.nbNodes(); ++i)
    {
        Integer aNodeId = m_mesh.nodeId(i);

        CMshNode<double> aNode = m_mesh.node(aNodeId);

        pnt[0] = aNode.x();
        pnt[1] = aNode.y();
        pnt[2] = aNode.z();

        Ng_AddPoint(mesh, pnt);
    }

    for (Integer i = 0; i < m_mesh.nbElements(); ++i)
    {
        Integer anElementId = m_mesh.elementId(i);

        CMshElement<double> anElement = m_mesh.element(anElementId);

        if (anElement.elementType() == ET_TRIANGLE)
        {
            tri[0] = anElement.nodeId(0) + 1;
            tri[1] = anElement.nodeId(1) + 1;
            tri[2] = anElement.nodeId(2) + 1;

            Ng_AddSurfaceElement(mesh, NG_TRIG, tri);
        }
    }

    QElapsedTimer timer;

    timer.restart();

    Ng_GenerateVolumeMesh(mesh, h);

    std::cout << "Elapsed time: " << timer.elapsed() / 1000.0 << " seconds" << std::endl;

    m_mesh.reset();

    int np = Ng_GetNP(mesh);
    for (int i = 1; i <= np; ++i)
    {
        Ng_GetPoint (mesh, i, pnt);

        CMshNode<double> aNode(pnt[0], pnt[1], pnt[2]);
        m_mesh.addNode(i - 1, aNode);
    }

    int ns = Ng_GetNSE(mesh);
    for (int i = 1; i <= ns; ++i)
    {
        Ng_GetSurfaceElement (mesh, i, tri);

        CMshElement<double> anElement(ET_TRIANGLE);
        anElement.addNodeId(tri[0] - 1);
        anElement.addNodeId(tri[1] - 1);
        anElement.addNodeId(tri[2] - 1);

        Integer anElementId = m_mesh.nextElementId();
        m_mesh.addElement(anElementId, anElement);
    }

    int ne = Ng_GetNE(mesh);
    for (int i = 1; i <= ne; ++i)
    {
        Ng_GetVolumeElement (mesh, i, tet);

        CMshElement<double> anElement(ET_TETRAHEDRON);
        anElement.addNodeId(tet[0] - 1);
        anElement.addNodeId(tet[1] - 1);
        anElement.addNodeId(tet[2] - 1);
        anElement.addNodeId(tet[3] - 1);

        Integer anElementId = m_mesh.nextElementId();
        m_mesh.addElement(anElementId, anElement);
    }

    double minq, maxq, aveq;

    m_mesh.meshQuality(minq, maxq, aveq);
    std::cout << "Final mesh quality: " << aveq << " [" << minq << ", " << maxq << "]" << std::endl;

    Ng_Exit();

    #else

    CMshTetrahedronMesher<double> aTetrahedronMesher;

    std::vector<CMshNode<double> > sInnerNodes;

    m_mesh.mergeNodes(tolerance);

    m_mesh.generateFaces(tolerance*tolerance);

    QElapsedTimer timer;

    timer.restart();

    aTetrahedronMesher.generate(m_mesh, 999999, sInnerNodes, h, 0.1, tolerance*tolerance);

    std::cout << "Elapsed time: " << timer.elapsed() / 1000.0 << " seconds" << std::endl;

    m_mesh.addMesh(aTetrahedronMesher.mesh());

    m_mesh.mergeNodes(tolerance);

    m_mesh.renumber();

    double minq, maxq, aveq;

    m_mesh.meshQuality(minq, maxq, aveq);
    std::cout << "Mesh quality: " << aveq << " [" << minq << ", " << maxq << "]" << std::endl;

    #endif

    /*
    #ifdef USE_NETGEN
    std::cout << "Performing optimization..." << std::endl;

    Ng_Init();

    Ng_Mesh* mesh = Ng_NewMesh();

    double pnt[3];

    for (Integer i = 0; i < m_mesh.nbNodes(); ++i)
    {
        Integer aNodeId = m_mesh.nodeId(i);

        CMshNode<double> aNode = m_mesh.node(aNodeId);

        pnt[0] = aNode.x();
        pnt[1] = aNode.y();
        pnt[2] = aNode.z();

        Ng_AddPoint(mesh, pnt);
    }

    int tri[3];

    for (Integer i = 0; i < m_mesh.nbElements(); ++i)
    {
        Integer anElementId = m_mesh.elementId(i);

        CMshElement<double> anElement = m_mesh.element(anElementId);

        if (anElement.elementType() == ET_TRIANGLE)
        {
            tri[0] = anElement.nodeId(0) + 1;
            tri[1] = anElement.nodeId(1) + 1;
            tri[2] = anElement.nodeId(2) + 1;

            Ng_AddSurfaceElement(mesh, NG_TRIG, tri);
        }
    }

    int tet[4];

    for (Integer i = 0; i < m_mesh.nbElements(); ++i)
    {
        Integer anElementId = m_mesh.elementId(i);

        CMshElement<double> anElement = m_mesh.element(anElementId);

        if (anElement.elementType() == ET_TETRAHEDRON)
        {
            tet[0] = anElement.nodeId(0) + 1;
            tet[1] = anElement.nodeId(1) + 1;
            tet[2] = anElement.nodeId(2) + 1;
            tet[3] = anElement.nodeId(3) + 1;

            Ng_AddVolumeElement(mesh, NG_TET, tet);
        }
    }

    Ng_OptimizeVolumeMesh(mesh, h);

    m_mesh.reset();

    int np = Ng_GetNP(mesh);
    for (int i = 1; i <= np; ++i)
    {
        Ng_GetPoint (mesh, i, pnt);

        CMshNode<double> aNode(pnt[0], pnt[1], pnt[2]);
        m_mesh.addNode(i - 1, aNode);
    }

    int ns = Ng_GetNSE(mesh);
    for (int i = 1; i <= ns; ++i)
    {
        Ng_GetSurfaceElement (mesh, i, tri);

        CMshElement<double> anElement(ET_TRIANGLE);
        anElement.addNodeId(tri[0] - 1);
        anElement.addNodeId(tri[1] - 1);
        anElement.addNodeId(tri[2] - 1);

        Integer anElementId = m_mesh.nextElementId();
        m_mesh.addElement(anElementId, anElement);
    }

    int ne = Ng_GetNE(mesh);
    for (int i = 1; i <= ne; ++i)
    {
        Ng_GetVolumeElement (mesh, i, tet);

        CMshElement<double> anElement(ET_TETRAHEDRON);
        anElement.addNodeId(tet[0] - 1);
        anElement.addNodeId(tet[1] - 1);
        anElement.addNodeId(tet[2] - 1);
        anElement.addNodeId(tet[3] - 1);

        Integer anElementId = m_mesh.nextElementId();
        m_mesh.addElement(anElementId, anElement);
    }

    m_mesh.meshQuality(minq, maxq, aveq);
    std::cout << "Final mesh quality: " << aveq << " [" << minq << ", " << maxq << "]" << std::endl;

    Ng_Exit();
    #endif
    */

    return true;

}

void OccMesher::stitchMesh()
{

    QElapsedTimer timer;

    timer.restart();

    m_mesh.generateFaces(1E-2);
    m_mesh.collapseNakedEdges(1E-2);
    m_mesh.mergeNodes(1E-2);
    m_mesh.generateFaces(1E-2);

    std::cout << "Elapsed time: " << timer.elapsed() / 1000.0 << " seconds" << std::endl;

}

