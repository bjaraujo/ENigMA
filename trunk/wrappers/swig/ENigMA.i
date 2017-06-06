%module(directors="1") ENigMA

// Remove some warnings
#pragma SWIG nowarn=362,503,401,389,516,511

// Use STL support
%include <std_vector.i>
%include <std_string.i>
%include <std_map.i>

#if SWIGPYTHON || SWIGRUBY
%include <std_complex.i>
#endif
// Use C99 int support
%include <stdint.i>

// Use exceptions
%include "exception.i"

// Customize exception handling
%exception {
  try {
    $action
  } catch( std::exception &ex ) {
    const size_t e_size = 10240;
    char error_msg[e_size];
// TODO this should be replaces with some try compile stuff

%#ifdef _MSC_VER
    _snprintf_s( error_msg, e_size, e_size, "Exception thrown in ENigMA $symname: %s", ex.what() );
%#else
    snprintf( error_msg, e_size, "Exception thrown in ENigMA $symname: %s", ex.what() );
%#endif

    SWIG_exception( SWIG_RuntimeError, error_msg );
  } catch( ... ) {
    SWIG_exception( SWIG_UnknownError, "Unknown exception thrown in ENigMA $symname" );
  }
}

%include "typemaps.i"

%apply int    *OUTPUT {int& nIterations}; 
%apply int    *OUTPUT {int& smallEdge}; 
%apply int    *OUTPUT {int& bigEdge}; 

%apply double *OUTPUT {double& aDistance}; 
%apply double *OUTPUT {double& aTolerance}; 
%apply double *OUTPUT {double& aValue}; 
%apply double *OUTPUT {double& aValue1}; 
%apply double *OUTPUT {double& aValue2}; 
%apply double *OUTPUT {double& aMassError}; 

%apply double *OUTPUT {double& aMinQ}; 
%apply double *OUTPUT {double& aMaxQ}; 
%apply double *OUTPUT {double& aAveQ}; 

%apply double *OUTPUT {double& volumeFractionAct}; 
%apply double *OUTPUT {double& T}; 
%apply double *OUTPUT {double& d}; 
%apply double *OUTPUT {double& ru}; 
%apply double *OUTPUT {double& rv}; 
%apply double *OUTPUT {double& rw}; 
%apply double *OUTPUT {double& rp}; 
%apply double *OUTPUT {double& rT}; 
%apply double *OUTPUT {double& rs}; 

%apply double *OUTPUT {double& dmin}; 
%apply double *OUTPUT {double& dmax}; 
%apply double *OUTPUT {double& q}; 

#if !defined(SWIGLUA) && !defined(SWIGR)
%rename(Plus) operator +;
%rename(Minus) operator -;
%rename(Multiply) operator *;
%rename(Divide) operator /;
#endif

#if SWIGCSHARP
namespace std
{
  %template(StdVectorInt) vector<int>;
  %template(StdVectorFloat) vector<float>;
  %template(StdVectorDouble) vector<double>; 
}
#endif

%{
#include "CmnTypes.hpp"
#include "AnaFunction.hpp"
#include "AnaTemperature.hpp"
#include "GeoCoordinateSystem.hpp"
#include "GeoCoordinate.hpp"
#include "GeoVector.hpp"
#include "GeoNormal.hpp"
#include "GeoPlane.hpp"
#include "GeoBoundingBox.hpp"
#include "GeoVertexList.hpp"
#include "GeoLine.hpp"
#include "GeoLineList.hpp"
#include "GeoPolyline.hpp"
#include "GeoTriangle.hpp"
#include "GeoQuadrilateral.hpp"
#include "GeoTetrahedron.hpp"
#include "GeoTriangularPrism.hpp"
#include "GeoHexahedron.hpp"
#include "GeoPolyhedron.hpp"
#include "GeoAdtree.hpp"
#include "GeoHashGrid.hpp"
#include "GeoOctree.hpp"
#include "GeoRtree.hpp"
#include "StlFile.hpp"
#include "StlUtils.hpp"
#include "MshNode.hpp"
#include "MshFace.hpp"
#include "MshElement.hpp"
#include "MshMesh.hpp"
#include "MshBasicMesher.hpp"
#include "MshExtrudedMesher.hpp"
#include "MshTriangleMesher.hpp"
#include "MshQuadrilateralMesher.hpp"
#include "SleSystem.hpp"
#include "MatMaterial.hpp"
#include "PdeField.hpp"
#include "PdeEquation.hpp"
#include "PdeBoundaryCondition.hpp"
#include "FemBeam.hpp"
#include "FemTriangle.hpp"
#include "FemQuadrilateral.hpp"
#include "FemTetrahedron.hpp"
#include "FemTriangularPrism.hpp"
#include "FemHexahedron.hpp"
#include "FvmNode.hpp"
#include "FvmFace.hpp"
#include "FvmCell.hpp"
#include "FvmControlVolume.hpp"
#include "FvmMesh.hpp"
#include "FvmPisoSolver.hpp"
#include "FvmTemperatureSolver.hpp"
#include "FvmVofSolver.hpp"
#include "SphKernel.hpp"
#include "SphConvex.hpp"
#include "SphCubicSpline.hpp"
#include "SphGaussian.hpp"
#include "SphQuintic.hpp"
#include "SphSpiky.hpp"
#include "SphParticles.hpp"
#include "PosGmsh.hpp"
#include "PosVtk.hpp"
%}

%template(StdVectorCGeoPolygonDouble) std::vector<ENigMA::geometry::CGeoPolygon<double> >;
%template(StdVectorCGeoCoordinateDouble) std::vector<ENigMA::geometry::CGeoCoordinate<double> >;
%template(StdVectorCGeoTriangleDouble) std::vector<ENigMA::geometry::CGeoTriangle<double> >;
%template(StdVectorCGeoTetrahedronDouble) std::vector<ENigMA::geometry::CGeoTetrahedron<double> >;
%template(StdVectorCMshFaceDouble) std::vector<ENigMA::mesh::CMshFace<double> >;

%include "CmnTypes.hpp"

// Coordinate
%include "GeoCoordinate.hpp"

%template(CGeoCoordinateDouble) ENigMA::geometry::CGeoCoordinate<double>;

%extend ENigMA::geometry::CGeoCoordinate<double> {

    double x() { 
        return (*$self).x();
    }

    double y() { 
        return (*$self).y();
    }

    double z() { 
        return (*$self).z();
    }

    void setX(const double x) { 
        (*$self).x() = x;
    }

    void setY(const double y) { 
        (*$self).y() = y;
    }

    void setZ(const double z) { 
        (*$self).z() = z;
    }

    ENigMA::geometry::CGeoCoordinate<double> operator+(const ENigMA::geometry::CGeoCoordinate<double>& c) const {
        return (*$self) + c;
    }

    ENigMA::geometry::CGeoCoordinate<double> operator-(const ENigMA::geometry::CGeoCoordinate<double>& c) const {
        return (*$self) - c;
    }

    ENigMA::geometry::CGeoCoordinate<double> operator*(const double c) const {
        return (*$self) * c;
    }

    ENigMA::geometry::CGeoCoordinate<double> operator/(const double c) const {
        return (*$self) / c;
    }

};

// Coordinate System
%include "GeoCoordinateSystem.hpp"

%template(CGeoCoordinateSystemDouble) ENigMA::geometry::CGeoCoordinateSystem<double>;

%extend ENigMA::geometry::CGeoCoordinateSystem<double> {

    double operator()(const int nRow, const int nCol) {
        return (*$self)(nRow, nCol);
    }

    void operator()(const int nRow, const int nCol, const double value) {
        (*$self)(nRow, nCol) = value;
    }

};

// Vector
%include "GeoVector.hpp"

%template(CGeoVectorDouble) ENigMA::geometry::CGeoVector<double>;

%extend ENigMA::geometry::CGeoVector<double> {

    double x() { 
        return (*$self).x();
    }

    double y() { 
        return (*$self).y();
    }

    double z() { 
        return (*$self).z();
    }

    ENigMA::geometry::CGeoVector<double> operator+(const ENigMA::geometry::CGeoVector<double>& c) const {
        return (*$self) + c;
    }

    ENigMA::geometry::CGeoVector<double> operator-(const ENigMA::geometry::CGeoVector<double>& c) const {
        return (*$self) - c;
    }

    ENigMA::geometry::CGeoVector<double> operator*(const double c) const {
        return (*$self) * c;
    }

    ENigMA::geometry::CGeoVector<double> operator/(const double c) const {
        return (*$self) / c;
    }

};

// Normal
%include "GeoNormal.hpp"

%template(CGeoNormalDouble) ENigMA::geometry::CGeoNormal<double>;

%extend ENigMA::geometry::CGeoNormal<double> {

    double x() { 
        return (*$self).x();
    }

    double y() { 
        return (*$self).y();
    }

    double z() { 
        return (*$self).z();
    }

    void normalize() { 
        (*$self).normalize();
    }

    ENigMA::geometry::CGeoNormal<double> operator+(const ENigMA::geometry::CGeoNormal<double>& c) const {
        return (*$self) + c;
    }

    ENigMA::geometry::CGeoNormal<double> operator-(const ENigMA::geometry::CGeoNormal<double>& c) const {
        return (*$self) - c;
    }

    ENigMA::geometry::CGeoNormal<double> operator*(const double c) const {
        return (*$self) * c;
    }

    ENigMA::geometry::CGeoNormal<double> operator/(const double c) const {
        return (*$self) / c;
    }

};

// Bounding Box
%include "GeoBoundingBox.hpp"

%template(CGeoBoundingBoxDouble) ENigMA::geometry::CGeoBoundingBox<double>;

// Plane
%include "GeoPlane.hpp"

%template(CGeoPlaneDouble) ENigMA::geometry::CGeoPlane<double>;

%extend ENigMA::geometry::CGeoPlane<double> {

};

// Vertex List
%include "GeoVertexList.hpp"

%template(CGeoVertexListDouble) ENigMA::geometry::CGeoVertexList<double>;

// Line
%include "GeoLine.hpp"

%template(CGeoLineDouble) ENigMA::geometry::CGeoLine<double>;

%extend ENigMA::geometry::CGeoLine<double> {
    
    double length() { 
        return (*$self).length();
    }

};

// Line List
%include "GeoLineList.hpp"

%template(CGeoLineListDouble) ENigMA::geometry::CGeoLineList<double>;

// Polyline
%include "GeoPolyline.hpp"

%template(CGeoPolylineDouble) ENigMA::geometry::CGeoPolyline<double>;

%extend ENigMA::geometry::CGeoPolyline<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

%extend ENigMA::geometry::CGeoPolyline<double> {
    
    int nbVertices() { 
        return (*$self).nbVertices();
    }

};

// Triangle
%include "GeoTriangle.hpp"

%template(CGeoTriangleDouble) ENigMA::geometry::CGeoTriangle<double>;

%extend ENigMA::geometry::CGeoTriangle<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

%extend ENigMA::geometry::CGeoTriangle<double> {
    
    int nbVertices() { 
        return (*$self).nbVertices();
    }

};

// Quadrilateral
%include "GeoQuadrilateral.hpp"

%template(CGeoQuadrilateralDouble) ENigMA::geometry::CGeoQuadrilateral<double>;

%extend ENigMA::geometry::CGeoQuadrilateral<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

%extend ENigMA::geometry::CGeoQuadrilateral<double> {
    
    int nbVertices() { 
        return (*$self).nbVertices();
    }

};

// Tetrahedron
%include "GeoTetrahedron.hpp"

%template(CGeoTetrahedronDouble) ENigMA::geometry::CGeoTetrahedron<double>;

%extend ENigMA::geometry::CGeoTetrahedron<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

%extend ENigMA::geometry::CGeoTetrahedron<double> {
    
    int nbVertices() { 
        return (*$self).nbVertices();
    }

};

// Triangular Prism
%include "GeoTriangularPrism.hpp"

%template(CGeoTriangularPrismDouble) ENigMA::geometry::CGeoTriangularPrism<double>;

%extend ENigMA::geometry::CGeoTriangularPrism<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

// Hexahedron
%include "GeoHexahedron.hpp"

%template(CGeoHexahedronDouble) ENigMA::geometry::CGeoHexahedron<double>;

%extend ENigMA::geometry::CGeoHexahedron<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

%extend ENigMA::geometry::CGeoHexahedron<double> {
    
    int nbVertices() { 
        return (*$self).nbVertices();
    }

};

// Polygon
%include "GeoPolygon.hpp"

%template(CGeoPolygonDouble) ENigMA::geometry::CGeoPolygon<double>;

// Polyhedron
%include "GeoPolyhedron.hpp"

%template(CGeoPolyhedronDouble) ENigMA::geometry::CGeoPolyhedron<double>;

// Ad-Tree
%include "GeoAdtree.hpp"

%template(CGeoAdtreeDouble) ENigMA::geometry::CGeoAdtree<double>;

// Hash Grid
%include "GeoHashGrid.hpp"

%template(CGeoHashGridDouble) ENigMA::geometry::CGeoHashGrid<double>;

// Octree
%include "GeoOctree.hpp"

%template(CGeoOctreeDouble) ENigMA::geometry::CGeoOctree<double>;

// R-Tree
%include "GeoRtree.hpp"

%template(CGeoRtreeDouble) ENigMA::geometry::CGeoRtree<double>;

// STL - File
%include "StlFile.hpp"

%ignore ENigMA::stl::CStlStats<double>::fileSize;

%template(CStlEdgeDouble) ENigMA::stl::CStlEdge<double>;
%template(CStlFacetDouble) ENigMA::stl::CStlFacet<double>;
%template(CStlStatsDouble) ENigMA::stl::CStlStats<double>;
%template(CStlFileDouble) ENigMA::stl::CStlFile<double>;

// STL - Utils
%include "StlUtils.hpp"

%template(CStlUtilsDouble) ENigMA::stl::CStlUtils<double>;

// Node
%include "MshNode.hpp"

%template(CMshNodeDouble) ENigMA::mesh::CMshNode<double>;

%extend ENigMA::mesh::CMshNode<double> {

    double x() { 
        return (*$self).x();
    }

    double y() { 
        return (*$self).y();
    }

    double z() { 
        return (*$self).z();
    }

};

// Face
%include "MshFace.hpp"

%template(CMshFaceDouble) ENigMA::mesh::CMshFace<double>;

// Element
%include "MshElement.hpp"

%template(CMshElementDouble) ENigMA::mesh::CMshElement<double>;

// Mesh
%include "MshMesh.hpp"

%template(CMshMeshDouble) ENigMA::mesh::CMshMesh<double>;

// Basic Mesher
%include "MshBasicMesher.hpp"

%template(CMshBasicMesherDouble) ENigMA::mesh::CMshBasicMesher<double>;

// Extruded Mesher
%include "MshExtrudedMesher.hpp"

%template(CMshExtrudedMesherDouble) ENigMA::mesh::CMshExtrudedMesher<double>;

// Triangle Mesher
%include "MshTriangleMesher.hpp"

%ignore ENigMA::mesh::CMshTriangleMesher<double>::onUpdate;

%template(CMshTriangleMesherDouble) ENigMA::mesh::CMshTriangleMesher<double>;

// Quadrilateral Mesher
%include "MshQuadrilateralMesher.hpp"

%ignore ENigMA::mesh::CMshQuadrilateralMesher<double>::onUpdate;

%template(CMshQuadrilateralMesherDouble) ENigMA::mesh::CMshQuadrilateralMesher<double>;

// System of Linear Equations
%include "SleSystem.hpp"

%ignore ENigMA::sle::CSleSystem<double>::matrixA;
%ignore ENigMA::sle::CSleSystem<double>::vectorB;
%ignore ENigMA::sle::CSleSystem<double>::solve();

%template(CSleSystemDouble) ENigMA::sle::CSleSystem<double>;

%extend ENigMA::sle::CSleSystem<double> {

    ENigMA::sle::CSleSystem<double>& operator==(const double c) {
        return (*$self) = c;
    }

    ENigMA::sle::CSleSystem<double> operator+(const ENigMA::sle::CSleSystem<double>& c) const {
        return (*$self) + c;
    }

    ENigMA::sle::CSleSystem<double> operator-(const ENigMA::sle::CSleSystem<double>& c) const {
        return (*$self) - c;
    }

    ENigMA::sle::CSleSystem<double> operator*(const double c) const {
        return c * (*$self);
    }

    ENigMA::sle::CSleSystem<double> operator/(const double c) const {
        return 1.0 / c * (*$self);
    }

    ENigMA::sle::CSleSystem<double> operator-(const CSleSystem<double>& left, const CSleSystem<double>& right) {
        return left - right;
    }

    ENigMA::sle::CSleSystem<double> operator*(const double c, const CSleSystem<double>& right) {
        return c * right;
    }

}

// Analytical function
%include "AnaFunction.hpp"

%template(CAnaFunctionDouble) ENigMA::analytical::CAnaFunction<double>;

// Analytical temperature
%include "AnaTemperature.hpp"

%template(CAnaTemperatureDouble) ENigMA::analytical::CAnaTemperature<double>;

// Material
%include "MatMaterial.hpp"

%template(CMatMaterialDouble) ENigMA::material::CMatMaterial<double>;

// PDE Field
%include "PdeField.hpp"

%ignore ENigMA::pde::CPdeField<double>::u;
%ignore ENigMA::pde::CPdeField<double>::uFixed;
%ignore ENigMA::pde::CPdeField<double>::uSource;

%template(CPdeFieldDouble) ENigMA::pde::CPdeField<double>;

// PDE Equation
%include "PdeEquation.hpp"

%template(CPdeEquationDouble) ENigMA::pde::CPdeEquation<double>;

%template(ddt) ENigMA::pde::ddt<double>;
%template(laplacian) ENigMA::pde::laplacian<double>;
%template(divergence) ENigMA::pde::divergence<double>;

// Boundary condition
%include "PdeBoundaryCondition.hpp"

%template(CPdeBoundaryConditionDouble) ENigMA::pde::CPdeBoundaryCondition<double>;

/*
// FEM Beam
%include "FemBeam.hpp"

%template(CFemBeamDouble211) ENigMA::fem::CFemBeam<double, 2, 1, 1>;

// FEM Triangle
%include "FemTriangle.hpp"

%template(CFemTriangleDouble311) ENigMA::fem::CFemTriangle<double, 3, 1, 1>;

// FEM Quadrilateral
%include "FemQuadrilateral.hpp"

%template(CFemQuadrilateralDouble411) ENigMA::fem::CFemQuadrilateral<double, 4, 1, 1>;

// FEM Tetrahedron
%include "FemTetrahedron.hpp"

%template(CFemTetrahedronDouble411) ENigMA::fem::CFemTetrahedron<double, 4, 1, 1>;

// FEM Triangular Prism
%include "FemTriangularPrism.hpp"

%template(CFemTriangularPrismDouble611) ENigMA::fem::CFemTriangularPrism<double, 6, 1, 1>;

// FEM Hexahedron
%include "FemHexahedron.hpp"

%template(CFemHexahedronDouble811) ENigMA::fem::CFemHexahedron<double, 8, 1, 1>;
*/

// FVM Node
%include "FvmNode.hpp"

%template(CFvmNodeDouble) ENigMA::fvm::CFvmNode<double>;

// Fvm Face
%include "FvmFace.hpp"

%template(CFvmFaceDouble) ENigMA::fvm::CFvmFace<double>;

%extend ENigMA::fvm::CFvmFace<double> {
    
    ENigMA::geometry::CGeoCoordinate<double>& centroid() { 
        return (*$self).centroid();
    }

    void calculateCentroid() { 
        (*$self).calculateCentroid();
    }

};

// Fvm Cell
%include "FvmCell.hpp"

%template(CFvmCellDouble) ENigMA::fvm::CFvmCell<double>;

// Fvm Control Volume
%include "FvmControlVolume.hpp"

%template(CFvmControlVolumDouble) ENigMA::fvm::CFvmControlVolume<double>;

// FVM Mesh
%include "FvmMesh.hpp"

%template(CFvmMeshDouble) ENigMA::fvm::CFvmMesh<double>;

// FVM Piso Solver
%include "FvmPisoSolver.hpp"

%template(CFvmPisoSolverDouble) ENigMA::fvm::CFvmPisoSolver<double>;

// FVM Temperature Solver
%include "FvmTemperatureSolver.hpp"

%template(CFvmTemperatureSolverDouble) ENigMA::fvm::CFvmTemperatureSolver<double>;

// FVM Vof Solver
%include "FvmVofSolver.hpp"

%template(CFvmVofSolverDouble) ENigMA::fvm::CFvmVofSolver<double>;

// SPH Kernel
%include "SphKernel.hpp"

%template(CSphKernelDouble) ENigMA::sph::CSphKernel<double>;

// SPH Convex
%include "SphConvex.hpp"

%template(CSphConvexDouble) ENigMA::sph::CSphConvex<double>;

// SPH Cubic Spline
%include "SphCubicSpline.hpp"

%template(CSphCubicSplineDouble) ENigMA::sph::CSphCubicSpline<double>;

// SPH Gaussian
%include "SphGaussian.hpp"

%template(CSphGaussianDouble) ENigMA::sph::CSphGaussian<double>;

// SPH Quintic
%include "SphQuintic.hpp"

%template(CSphQuinticDouble) ENigMA::sph::CSphQuintic<double>;

// SPH Spiky
%include "SphSpiky.hpp"

%template(CSphSpikyDouble) ENigMA::sph::CSphSpiky<double>;

// SPH 

%include "SphParticles.hpp"

%template(CSphParticlesDouble) ENigMA::sph::CSphParticles<double>;

// Gmsh
%include "PosGmsh.hpp"

%template(CPosGmshDouble) ENigMA::post::CPosGmsh<double>;

// Vtk
%include "PosVtk.hpp"

%template(CPosVtkDouble) ENigMA::post::CPosVtk<double>;
