
project(ENigMA)

set(GEOMETRY_HEADERS
geometry/GeoCentroid.hpp
geometry/GeoCentroid_Imp.hpp
geometry/GeoArea.hpp
geometry/GeoArea_Imp.hpp
geometry/GeoCoordinate.hpp
geometry/GeoCoordinate_Imp.hpp
geometry/GeoBoundingBox.hpp
geometry/GeoBoundingBox_Imp.hpp
geometry/GeoCoordinateSystem.hpp
geometry/GeoHexahedron.hpp
geometry/GeoHexahedron_Imp.hpp
geometry/GeoLength.hpp
geometry/GeoLength_Imp.hpp
geometry/GeoLine.hpp
geometry/GeoLine_Imp.hpp
geometry/GeoNormal.hpp
geometry/GeoPlane.hpp
geometry/GeoPlane_Imp.hpp
geometry/GeoPolygon.hpp
geometry/GeoPolygon_Imp.hpp
geometry/GeoPolyline.hpp
geometry/GeoPolyline_Imp.hpp
geometry/GeoPolyhedron.hpp
geometry/GeoPolyhedron_Imp.hpp
geometry/GeoQuadrilateral.hpp
geometry/GeoQuadrilateral_Imp.hpp
geometry/GeoTetrahedron.hpp
geometry/GeoTetrahedron_Imp.hpp
geometry/GeoTriangle.hpp
geometry/GeoTriangle_Imp.hpp
geometry/GeoTriangularPrism.hpp
geometry/GeoTriangularPrism_Imp.hpp
geometry/GeoCircle.hpp
geometry/GeoCircle_Imp.hpp
geometry/GeoSphere.hpp
geometry/GeoSphere_Imp.hpp
geometry/GeoVector.hpp
geometry/GeoVector_Imp.hpp
geometry/GeoVertexList.hpp
geometry/GeoVertexList_Imp.hpp
geometry/GeoLineList.hpp
geometry/GeoLineList_Imp.hpp
geometry/GeoVolume.hpp
geometry/GeoVolume_Imp.hpp
geometry/GeoContainer.hpp
geometry/GeoContainer_Imp.hpp
geometry/GeoHashGrid.hpp
geometry/GeoHashGrid_Imp.hpp
geometry/GeoOctreeNode.hpp
geometry/GeoOctreeNode_Imp.hpp
geometry/GeoOctree.hpp
geometry/GeoOctree_Imp.hpp
geometry/GeoRtree.hpp
geometry/GeoRtree_Imp.hpp
geometry/GeoAdtree.hpp
geometry/GeoAdtree_Imp.hpp
geometry/GeoConvexHull.hpp
geometry/GeoConvexHull_Imp.hpp
geometry/GeoIntersectionType.hpp
geometry/GeoIntersectionType_Imp.hpp
geometry/GeoPredicates.hpp
geometry/GeoPredicates_Imp.hpp
)

set(MESH_HEADERS
mesh/MshNode.hpp
mesh/MshFace.hpp
mesh/MshFace_Imp.hpp
mesh/MshElement.hpp
mesh/MshElement_Imp.hpp
mesh/MshMesh.hpp
mesh/MshMesh_Imp.hpp
mesh/MshCheckMesh.hpp
mesh/MshCheckMesh_Imp.hpp
mesh/MshMeshQuery.hpp
mesh/MshMeshQuery_Imp.hpp
mesh/MshBasicMesher.hpp
mesh/MshBasicMesher_Imp.hpp
mesh/MshExtrudedMesher.hpp
mesh/MshExtrudedMesher_Imp.hpp
mesh/MshTriangle.hpp
mesh/MshTriangle_Imp.hpp
mesh/MshQuadrilateral.hpp
mesh/MshQuadrilateral_Imp.hpp
mesh/MshTetrahedron.hpp
mesh/MshTetrahedron_Imp.hpp
mesh/MshTriangularPrism.hpp
mesh/MshTriangularPrism_Imp.hpp
mesh/MshTriangleMesher.hpp
mesh/MshTriangleMesher_Imp.hpp
mesh/MshQuadrilateralMesher.hpp
mesh/MshQuadrilateralMesher_Imp.hpp
mesh/MshTetrahedronMesher.hpp
mesh/MshTetrahedronMesher_Imp.hpp
)

set(INT_HEADERS
integration/IntGaussIntegration.hpp
integration/IntHexahedron.hpp
integration/IntHexahedron_Imp.hpp
integration/IntQuadrilateral.hpp
integration/IntQuadrilateral_Imp.hpp
integration/IntTetrahedron.hpp
integration/IntTetrahedron_Imp.hpp
integration/IntTriangle.hpp
integration/IntTriangle_Imp.hpp
integration/IntTriangularPrism.hpp
integration/IntTriangularPrism_Imp.hpp
)

set(BEM_HEADERS
bem/BemTriangle.hpp
bem/BemTriangle_Imp.hpp
)

set(FDM_HEADERS
fdm/FdmGrid.hpp
fdm/FdmGrid_Imp.hpp
)

set(FEM_HEADERS
fem/FemElement.hpp
fem/FemElement_Imp.hpp
fem/FemEdge.hpp
fem/FemEdge_Imp.hpp
fem/FemFace.hpp
fem/FemFace_Imp.hpp
fem/FemSolid.hpp
fem/FemBeam.hpp
fem/FemBeam_Imp.hpp
fem/FemTriangle.hpp
fem/FemTriangle_Imp.hpp
fem/FemQuadrilateral.hpp
fem/FemQuadrilateral_Imp.hpp
fem/FemTetrahedron.hpp
fem/FemTetrahedron_Imp.hpp
fem/FemHexahedron.hpp
fem/FemHexahedron_Imp.hpp
fem/FemTriangularPrism.hpp
fem/FemTriangularPrism_Imp.hpp
)

set(CVFEM_HEADERS
cvfem/CvFemMesh.hpp
cvfem/CvFemMesh_Imp.hpp
)

set(FEM_THERMAL_HEADERS
fem/thermal/FemThermalElement.hpp
fem/thermal/FemThermalElement_Imp.hpp
fem/thermal/FemLinearTemperatureBeam.hpp
fem/thermal/FemLinearTemperatureBeam_Imp.hpp
fem/thermal/FemLinearTemperatureTriangle.hpp
fem/thermal/FemLinearTemperatureTriangle_Imp.hpp
fem/thermal/FemLinearTemperatureQuadrilateral.hpp
fem/thermal/FemLinearTemperatureQuadrilateral_Imp.hpp
fem/thermal/FemLinearTemperatureTetrahedron.hpp
fem/thermal/FemLinearTemperatureTetrahedron_Imp.hpp
fem/thermal/FemLinearTemperatureTriangularPrism.hpp
fem/thermal/FemLinearTemperatureTriangularPrism_Imp.hpp
fem/thermal/FemLinearTemperatureHexahedron.hpp
fem/thermal/FemLinearTemperatureHexahedron_Imp.hpp
)

set(FEM_STRUCTURAL_HEADERS
fem/structural/FemStructuralElement.hpp
fem/structural/FemStructuralElement_Imp.hpp
fem/structural/FemConstantStrainTriangle.hpp
fem/structural/FemConstantStrainTriangle_Imp.hpp
fem/structural/FemConstantStrainTetrahedron.hpp
fem/structural/FemConstantStrainTetrahedron_Imp.hpp
)

set(FEM_FLOW_HEADERS
fem/flow/FemFlowElement.hpp
fem/flow/FemFlowElement_Imp.hpp
fem/flow/FemFlowTriangle.hpp
fem/flow/FemFlowTriangle_Imp.hpp
fem/flow/FemFlowTetrahedron.hpp
fem/flow/FemFlowTetrahedron_Imp.hpp
)

set(FVM_HEADERS
fvm/FvmNode.hpp
fvm/FvmFace.hpp
fvm/FvmFace_Imp.hpp
fvm/FvmCell.hpp
fvm/FvmCell_Imp.hpp
fvm/FvmControlVolume.hpp
fvm/FvmControlVolume_Imp.hpp
fvm/FvmMesh.hpp
fvm/FvmMesh_Imp.hpp
fvm/FvmMeshSearch.hpp
fvm/FvmMeshSearch_Imp.hpp
)

set(FVM_THERMAL_HEADERS
fvm/FvmTemperatureSolver.hpp
fvm/FvmTemperatureSolver_Imp.hpp
)

set(SPH_HEADERS
sph/SphKernel.hpp
sph/SphKernel_Imp.hpp
sph/SphCubicSpline.hpp
sph/SphCubicSpline_Imp.hpp
sph/SphGaussian.hpp
sph/SphGaussian_Imp.hpp
sph/SphQuintic.hpp
sph/SphQuintic_Imp.hpp
sph/SphSpiky.hpp
sph/SphSpiky_Imp.hpp
sph/SphConvex.hpp
sph/SphConvex_Imp.hpp
)

set(SLE_HEADERS
sle/SleSystem.hpp
sle/SleSystem_Imp.hpp
)

set(PDE_HEADERS
pde/PdeEquation.hpp
pde/PdeEquation_Imp.hpp
pde/PdeField.hpp
pde/PdeField_Imp.hpp
pde/PdeBoundaryCondition.hpp
pde/PdeBoundaryCondition_Imp.hpp
)

set(PDE_FDM_HEADERS
pde/fdm/FdmOperators.hpp
pde/fdm/FdmOperators_Imp.hpp
)

set(PDE_FDM_GENERIC_HEADERS
pde/fdm/generic/FdmGenericOperators_Imp.hpp
)

set(PDE_FDM_THERMAL_HEADERS
pde/fdm/thermal/FdmThermalOperators_Imp.hpp
)

set(PDE_BEM_HEADERS
pde/bem/BemOperators.hpp
pde/bem/BemOperators_Imp.hpp
)

set(PDE_BEM_THERMAL_HEADERS
pde/bem/thermal/BemThermalOperators_Imp.hpp
)

set(PDE_FEM_HEADERS
pde/fem/FemOperators.hpp
pde/fem/FemOperators_Imp.hpp
)

set(PDE_FEM_GENERIC_HEADERS
pde/fem/generic/FemGenericOperators_Imp.hpp
)

set(PDE_FEM_THERMAL_HEADERS
pde/fem/thermal/FemThermalOperators_Imp.hpp
)

set(PDE_FEM_STRUCTURAL_HEADERS
pde/fem/structural/FemStructuralOperators_Imp.hpp
)

set(PDE_FEM_FLOW_HEADERS
pde/fem/flow/FemFlowOperators_Imp.hpp
)

set(PDE_FVM_HEADERS
pde/fvm/FvmOperators.hpp
pde/fvm/FvmOperators_Imp.hpp
)

set(PDE_FVM_GENERIC_HEADERS
pde/fvm/generic/FvmGenericOperators_Imp.hpp
)

set(PDE_FVM_THERMAL_HEADERS
pde/fvm/thermal/FvmThermalOperators_Imp.hpp
)

set(PDE_FVM_FLOW_HEADERS
pde/fvm/flow/FvmFlowOperators_Imp.hpp
)

set(MAT_HEADERS
material/MatMaterial.hpp
material/MatMaterial_Imp.hpp
)

set(ANA_HEADERS
analytical/AnaFunction.hpp
analytical/AnaFunction_Imp.hpp
analytical/AnaTemperature.hpp
analytical/AnaTemperature_Imp.hpp
)

set(POS_HEADERS
post/PosGmsh.hpp
post/PosGmsh_Imp.hpp
post/PosVtk.hpp
post/PosVtk_Imp.hpp
post/PosGnuplot.hpp
post/PosGnuplot_Imp.hpp
post/PosQuickMesh.hpp
post/PosQuickMesh_Imp.hpp
)

set(STL_HEADERS
stl/StlUtils.hpp
stl/StlUtils_Imp.hpp
stl/StlFile.hpp
stl/StlFile_Imp.hpp
)

set(SOL_HEADERS
solvers/FemCbsSolver.hpp
solvers/FemCbsSolver_Imp.hpp
solvers/FvmPisoSolver.hpp
solvers/FvmPisoSolver_Imp.hpp
solvers/FvmVofSolver.hpp
solvers/FvmVofSolver_Imp.hpp
solvers/FvmTemperatureSolver.hpp
solvers/FvmTemperatureSolver_Imp.hpp
solvers/LbmLidDrivenSolver.hpp
solvers/LbmLidDrivenSolver_Imp.hpp
solvers/SphParticles.hpp
solvers/SphParticles_Imp.hpp
)

source_group("Header Files\\geometry" FILES ${GEOMETRY_HEADERS})
source_group("Header Files\\mesh" FILES ${MESH_HEADERS})
source_group("Header Files\\material" FILES ${MAT_HEADERS})
source_group("Header Files\\bem" FILES ${BEM_HEADERS})
source_group("Header Files\\fdm" FILES ${FDM_HEADERS})
source_group("Header Files\\fem" FILES ${FEM_HEADERS})
source_group("Header Files\\fem\\thermal" FILES ${FEM_THERMAL_HEADERS})
source_group("Header Files\\fem\\structural" FILES ${FEM_STRUCTURAL_HEADERS})
source_group("Header Files\\fem\\flow" FILES ${FEM_FLOW_HEADERS})
source_group("Header Files\\fvm" FILES ${FVM_HEADERS})
source_group("Header Files\\fvm\\thermal" FILES ${FVM_THERMAL_HEADERS})
source_group("Header Files\\fvm\\flow" FILES ${FVM_FLOW_HEADERS})
source_group("Header Files\\sph" FILES ${SPH_HEADERS})
source_group("Header Files\\sle" FILES ${SLE_HEADERS})
source_group("Header Files\\pde" FILES ${PDE_HEADERS})
source_group("Header Files\\integration" FILES ${INT_HEADERS})
source_group("Header Files\\pde\\bem" FILES ${PDE_BEM_HEADERS})
source_group("Header Files\\pde\\bem\\thermal" FILES ${PDE_BEM_THERMAL_HEADERS})
source_group("Header Files\\pde\\fdm" FILES ${PDE_FDM_HEADERS})
source_group("Header Files\\pde\\fdm\\generic" FILES ${PDE_FDM_GENERIC_HEADERS})
source_group("Header Files\\pde\\fdm\\thermal" FILES ${PDE_FDM_THERMAL_HEADERS})
source_group("Header Files\\pde\\fem" FILES ${PDE_FEM_HEADERS})
source_group("Header Files\\pde\\fem\\generic" FILES ${PDE_FEM_GENERIC_HEADERS})
source_group("Header Files\\pde\\fem\\thermal" FILES ${PDE_FEM_THERMAL_HEADERS})
source_group("Header Files\\pde\\fem\\structural" FILES ${PDE_FEM_STRUCTURAL_HEADERS})
source_group("Header Files\\pde\\fem\\flow" FILES ${PDE_FEM_FLOW_HEADERS})
source_group("Header Files\\pde\\fvm" FILES ${PDE_FVM_HEADERS})
source_group("Header Files\\pde\\fvm\\generic" FILES ${PDE_FVM_GENERIC_HEADERS})
source_group("Header Files\\pde\\fvm\\thermal" FILES ${PDE_FVM_THERMAL_HEADERS})
source_group("Header Files\\analytical" FILES ${ANA_HEADERS})
source_group("Header Files\\post" FILES ${POS_HEADERS})
source_group("Header Files\\stl" FILES ${STL_HEADERS})
