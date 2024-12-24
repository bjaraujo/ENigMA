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

%{
#include "Eigen/Dense"
#include "Eigen/Sparse"
%}

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

#if !defined(SWIGLUA) && !defined(SWIGR) && !defined(SWIGPYTHON)
%rename(Equal) operator =;
%rename(PlusEqual) operator +=;
%rename(MinusEqual) operator -=;
%rename(MultiplyEqual) operator *=;
%rename(DivideEqual) operator /=;
%rename(Plus) operator +;
%rename(Minus) operator -;
%rename(Multiply) operator *;
%rename(Divide) operator /;
%rename(IndexIntoConst) operator[](unsigned idx) const;
%rename(IndexInto) operator[](unsigned idx);
%rename(Functor) operator ();
#endif

#if SWIGPYTHON
%rename(__add__) operator +;
%rename(__sub__) operator -;
%rename(__mul__) operator *;
%rename(__div__) operator /;
#endif

namespace std
{
  %template(StdVectorInt) vector<int>;
  %template(StdVectorFloat) vector<float>;
  %template(StdVectorDouble) vector<double>;
}

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
#include "MshTetrahedronMesher.hpp"
#include "SleSystem.hpp"
#include "MatMaterial.hpp"
#include "PdeField.hpp"
#include "PdeEquation.hpp"
#include "PdeBoundaryCondition.hpp"
#include "IntGaussIntegration.hpp"
#include "FemBeam.hpp"
#include "FemTriangle.hpp"
#include "FemQuadrilateral.hpp"
#include "FemTetrahedron.hpp"
#include "FemTriangularPrism.hpp"
#include "FemHexahedron.hpp"
#include "FemCbsSolver.hpp"
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

%template(StdVectorCGeoPolygon) std::vector<ENigMA::geometry::CGeoPolygon<double> >;
%template(StdVectorCGeoCoordinate) std::vector<ENigMA::geometry::CGeoCoordinate<double> >;
%template(StdVectorCGeoTriangle) std::vector<ENigMA::geometry::CGeoTriangle<double> >;
%template(StdVectorCGeoTetrahedron) std::vector<ENigMA::geometry::CGeoTetrahedron<double> >;
%template(StdVectorCMshFace) std::vector<ENigMA::mesh::CMshFace<double> >;

%include "CmnTypes.hpp"

// Coordinate
%include "GeoCoordinate.hpp"

%template(CGeoCoordinate) ENigMA::geometry::CGeoCoordinate<double>;

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

%template(CGeoCoordinateSystem) ENigMA::geometry::CGeoCoordinateSystem<double>;

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

%template(CGeoVector) ENigMA::geometry::CGeoVector<double>;

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

%template(CGeoNormal) ENigMA::geometry::CGeoNormal<double>;

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

%template(CGeoBoundingBox) ENigMA::geometry::CGeoBoundingBox<double>;

// Plane
%include "GeoPlane.hpp"

%template(CGeoPlane) ENigMA::geometry::CGeoPlane<double>;

%extend ENigMA::geometry::CGeoPlane<double> {

};

// Vertex List
%include "GeoVertexList.hpp"

%template(CGeoVertexList) ENigMA::geometry::CGeoVertexList<double>;

// Line
%include "GeoLine.hpp"

%template(CGeoLine) ENigMA::geometry::CGeoLine<double>;

%extend ENigMA::geometry::CGeoLine<double> {
    
    double length() { 
        return (*$self).length();
    }

};

// Line List
%include "GeoLineList.hpp"

%template(CGeoLineList) ENigMA::geometry::CGeoLineList<double>;

// Polyline
%include "GeoPolyline.hpp"

%template(CGeoPolyline) ENigMA::geometry::CGeoPolyline<double>;

%extend ENigMA::geometry::CGeoPolyline<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

// Triangle
%include "GeoTriangle.hpp"

%template(CGeoTriangle) ENigMA::geometry::CGeoTriangle<double>;

%extend ENigMA::geometry::CGeoTriangle<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

// Quadrilateral
%include "GeoQuadrilateral.hpp"

%template(CGeoQuadrilateral) ENigMA::geometry::CGeoQuadrilateral<double>;

%extend ENigMA::geometry::CGeoQuadrilateral<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

// Tetrahedron
%include "GeoTetrahedron.hpp"

%template(CGeoTetrahedron) ENigMA::geometry::CGeoTetrahedron<double>;

%extend ENigMA::geometry::CGeoTetrahedron<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

// Triangular Prism
%include "GeoTriangularPrism.hpp"

%template(CGeoTriangularPrism) ENigMA::geometry::CGeoTriangularPrism<double>;

%extend ENigMA::geometry::CGeoTriangularPrism<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

// Hexahedron
%include "GeoHexahedron.hpp"

%template(CGeoHexahedron) ENigMA::geometry::CGeoHexahedron<double>;

%extend ENigMA::geometry::CGeoHexahedron<double> {
    
    void addVertex(ENigMA::geometry::CGeoCoordinate<double>& vertex) { 
        (*$self).addVertex(vertex);
    }

};

// Polygon
%include "GeoPolygon.hpp"

%template(CGeoPolygon) ENigMA::geometry::CGeoPolygon<double>;

// Polyhedron
%include "GeoPolyhedron.hpp"

%template(CGeoPolyhedron) ENigMA::geometry::CGeoPolyhedron<double>;

// Ad-Tree
%include "GeoAdtree.hpp"

%template(CGeoAdtree) ENigMA::geometry::CGeoAdtree<double>;

// Hash Grid
%include "GeoHashGrid.hpp"

%template(CGeoHashGrid) ENigMA::geometry::CGeoHashGrid<double>;

%extend ENigMA::geometry::CGeoHashGrid<double> {
    
    void addGeometricObject(int id, ENigMA::geometry::CGeoCoordinate<double>& coordinate) { 
        (*$self).addGeometricObject(id, coordinate);
    }

};

// Octree
%include "GeoOctree.hpp"

%template(CGeoOctree) ENigMA::geometry::CGeoOctree<double>;

%extend ENigMA::geometry::CGeoOctree<double> {
    
    void addGeometricObject(int id, ENigMA::geometry::CGeoCoordinate<double>& coordinate) { 
        (*$self).addGeometricObject(id, coordinate);
    }

};

// R-Tree
%include "GeoRtree.hpp"

%template(CGeoRtree) ENigMA::geometry::CGeoRtree<double>;

// STL - File
%include "StlFile.hpp"

%ignore ENigMA::stl::CStlStats<double>::fileSize;

%template(CStlEdge) ENigMA::stl::CStlEdge<double>;
%template(CStlFacet) ENigMA::stl::CStlFacet<double>;
%template(CStlStats) ENigMA::stl::CStlStats<double>;
%template(CStlFile) ENigMA::stl::CStlFile<double>;

// STL - Utils
%include "StlUtils.hpp"

%template(CStlUtils) ENigMA::stl::CStlUtils<double>;

// Node
%include "MshNode.hpp"

%template(CMshNode) ENigMA::mesh::CMshNode<double>;

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

%template(CMshFace) ENigMA::mesh::CMshFace<double>;

// Element
%include "MshElement.hpp"

%template(CMshElement) ENigMA::mesh::CMshElement<double>;

// Mesh
%include "MshMesh.hpp"

%template(CMshMesh) ENigMA::mesh::CMshMesh<double>;

// Basic Mesher
%include "MshBasicMesher.hpp"

%template(CMshBasicMesher) ENigMA::mesh::CMshBasicMesher<double>;

// Extruded Mesher
%include "MshExtrudedMesher.hpp"

%template(CMshExtrudedMesher) ENigMA::mesh::CMshExtrudedMesher<double>;

// Triangle Mesher
%include "MshTriangleMesher.hpp"

%ignore ENigMA::mesh::CMshTriangleMesher<double>::onUpdate;

%template(CMshTriangleMesher) ENigMA::mesh::CMshTriangleMesher<double>;

// Quadrilateral Mesher
%include "MshQuadrilateralMesher.hpp"

%ignore ENigMA::mesh::CMshQuadrilateralMesher<double>::onUpdate;

%template(CMshQuadrilateralMesher) ENigMA::mesh::CMshQuadrilateralMesher<double>;

// Quadrilateral Mesher
%include "MshTetrahedronMesher.hpp"

%ignore ENigMA::mesh::CMshTetrahedronMesher<double>::onUpdate;

%template(CMshTetrahedronMesher) ENigMA::mesh::CMshTetrahedronMesher<double>;

// System of Linear Equations
%include "SleSystem.hpp"

%template(CSleSystem) ENigMA::sle::CSleSystem<double>;

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

    ENigMA::sle::CSleSystem<double> operator*(const Eigen::Matrix<double, Eigen::Dynamic, 1>& left, const CSleSystem<double>& right) {
        return left * right;
    }

    ENigMA::sle::CSleSystem<double>& setRhs(const double c) {
        return (*$self) = c;
    }
    
    ENigMA::sle::CSleSystem<double>& setRhs(Eigen::Matrix<double, Eigen::Dynamic, 1> right) {
        return (*$self) = right;
    }

}

// Analytical function
%include "AnaFunction.hpp"

%template(CAnaFunction) ENigMA::analytical::CAnaFunction<double>;

// Analytical temperature
%include "AnaTemperature.hpp"

%template(CAnaTemperature) ENigMA::analytical::CAnaTemperature<double>;

// Material
%include "MatMaterial.hpp"

%template(CMatMaterial) ENigMA::material::CMatMaterial<double>;

// PDE Field
%include "PdeField.hpp"

%ignore ENigMA::pde::CPdeField<double>::uFixed;
%ignore ENigMA::pde::CPdeField<double>::uSource;

%template(CPdeField) ENigMA::pde::CPdeField<double>;

%extend ENigMA::pde::CPdeField<double> {

    ENigMA::sle::CSleSystem<double> operator*(const ENigMA::sle::CSleSystem<double>& c) const {
        return (*$self).u * c;
    }

}

// PDE Equation
%include "PdeEquation.hpp"

%template(CPdeEquation) ENigMA::pde::CPdeEquation<double>;

%template(ddt) ENigMA::pde::ddt<double>;
%template(laplacian) ENigMA::pde::laplacian<double>;
%template(divergence) ENigMA::pde::divergence<double>;
%template(gradient) ENigMA::pde::gradient<double>;

// Boundary condition
%include "PdeBoundaryCondition.hpp"

%template(CPdeBoundaryCondition) ENigMA::pde::CPdeBoundaryCondition<double>;

/*
// FEM Beam
%include "FemBeam.hpp"

%template(CFemBeam) ENigMA::fem::CFemBeam<double, 2, 1, 1>;

// FEM Triangle
%include "FemTriangle.hpp"

%template(CFemTriangle) ENigMA::fem::CFemTriangle<double, 3, 1, 1>;

// FEM Quadrilateral
%include "FemQuadrilateral.hpp"

%template(CFemQuadrilateral) ENigMA::fem::CFemQuadrilateral<double, 4, 1, 1>;

// FEM Tetrahedron
%include "FemTetrahedron.hpp"

%template(CFemTetrahedron) ENigMA::fem::CFemTetrahedron<double, 4, 1, 1>;

// FEM Triangular Prism
%include "FemTriangularPrism.hpp"

%template(CFemTriangularPrism) ENigMA::fem::CFemTriangularPrism<double, 6, 1, 1>;

// FEM Hexahedron
%include "FemHexahedron.hpp"

%template(CFemHexahedron) ENigMA::fem::CFemHexahedron<double, 8, 1, 1>;
*/

// FEM CBS Solver
%include "FemCbsSolver.hpp"

%template(CFemCbsSolver2) ENigMA::fem::CFemCbsSolver<double, 2>;
%template(CFemCbsSolver3) ENigMA::fem::CFemCbsSolver<double, 3>;

// FVM Node
%include "FvmNode.hpp"

%template(CFvmNode) ENigMA::fvm::CFvmNode<double>;

// Fvm Face
%include "FvmFace.hpp"

%template(CFvmFace) ENigMA::fvm::CFvmFace<double>;

%extend ENigMA::fvm::CFvmFace<double> {
    
    ENigMA::geometry::CGeoCoordinate<double> centroid() { 
        return (*$self).centroid();
    }

    double area() { 
        return (*$self).area();
    }
};

// Fvm Cell
%include "FvmCell.hpp"

%template(CFvmCell) ENigMA::fvm::CFvmCell<double>;

// Fvm Control Volume
%include "FvmControlVolume.hpp"

%template(CFvmControlVolume) ENigMA::fvm::CFvmControlVolume<double>;

// FVM Mesh
%include "FvmMesh.hpp"

%template(CFvmMesh) ENigMA::fvm::CFvmMesh<double>;

// FVM Piso Solver
%include "FvmPisoSolver.hpp"

%template(CFvmPisoSolver) ENigMA::fvm::CFvmPisoSolver<double>;

// FVM Temperature Solver
%include "FvmTemperatureSolver.hpp"

%template(CFvmTemperatureSolver) ENigMA::fvm::CFvmTemperatureSolver<double>;

// FVM Vof Solver
%include "FvmVofSolver.hpp"

%template(CFvmVofSolver) ENigMA::fvm::CFvmVofSolver<double>;

// SPH Kernel
%include "SphKernel.hpp"

%template(CSphKernel) ENigMA::sph::CSphKernel<double>;

// SPH Convex
%include "SphConvex.hpp"

%template(CSphConvex) ENigMA::sph::CSphConvex<double>;

// SPH Cubic Spline
%include "SphCubicSpline.hpp"

%template(CSphCubicSpline) ENigMA::sph::CSphCubicSpline<double>;

// SPH Gaussian
%include "SphGaussian.hpp"

%template(CSphGaussian) ENigMA::sph::CSphGaussian<double>;

// SPH Quintic
%include "SphQuintic.hpp"

%template(CSphQuintic) ENigMA::sph::CSphQuintic<double>;

// SPH Spiky
%include "SphSpiky.hpp"

%template(CSphSpiky) ENigMA::sph::CSphSpiky<double>;

// SPH

%include "SphParticles.hpp"

%template(CSphParticles) ENigMA::sph::CSphParticles<double>;

// Gmsh
%include "PosGmsh.hpp"

%template(CPosGmsh) ENigMA::post::CPosGmsh<double>;

// Vtk
%include "PosVtk.hpp"

%template(CPosVtk) ENigMA::post::CPosVtk<double>;
