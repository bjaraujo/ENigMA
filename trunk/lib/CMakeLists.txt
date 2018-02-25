
project(ENigMALib CXX)

set(GEOMETRY_HEADERS
../src/geometry/GeoCentroid.hpp
../src/geometry/GeoCentroid_Imp.hpp
../src/geometry/GeoArea.hpp
../src/geometry/GeoArea_Imp.hpp
../src/geometry/GeoCoordinate.hpp
../src/geometry/GeoCoordinate_Imp.hpp
../src/geometry/GeoBoundingBox.hpp
../src/geometry/GeoBoundingBox_Imp.hpp
../src/geometry/GeoCoordinateSystem.hpp
../src/geometry/GeoHexahedron.hpp
../src/geometry/GeoHexahedron_Imp.hpp
../src/geometry/GeoLength.hpp
../src/geometry/GeoLength_Imp.hpp
../src/geometry/GeoLine.hpp
../src/geometry/GeoLine_Imp.hpp
../src/geometry/GeoNormal.hpp
../src/geometry/GeoPlane.hpp
../src/geometry/GeoPlane_Imp.hpp
../src/geometry/GeoPolygon.hpp
../src/geometry/GeoPolygon_Imp.hpp
../src/geometry/GeoPolyline.hpp
../src/geometry/GeoPolyline_Imp.hpp
../src/geometry/GeoPolyhedron.hpp
../src/geometry/GeoPolyhedron_Imp.hpp
../src/geometry/GeoQuadrilateral.hpp
../src/geometry/GeoQuadrilateral_Imp.hpp
../src/geometry/GeoTetrahedron.hpp
../src/geometry/GeoTetrahedron_Imp.hpp
../src/geometry/GeoTriangle.hpp
../src/geometry/GeoTriangle_Imp.hpp
../src/geometry/GeoTriangularPrism.hpp
../src/geometry/GeoTriangularPrism_Imp.hpp
../src/geometry/GeoCircle.hpp
../src/geometry/GeoCircle_Imp.hpp
../src/geometry/GeoSphere.hpp
../src/geometry/GeoSphere_Imp.hpp
../src/geometry/GeoVector.hpp
../src/geometry/GeoVector_Imp.hpp
../src/geometry/GeoVertexList.hpp
../src/geometry/GeoVertexList_Imp.hpp
../src/geometry/GeoLineList.hpp
../src/geometry/GeoLineList_Imp.hpp
../src/geometry/GeoVolume.hpp
../src/geometry/GeoVolume_Imp.hpp
../src/geometry/GeoContainer.hpp
../src/geometry/GeoContainer_Imp.hpp
../src/geometry/GeoHashGrid.hpp
../src/geometry/GeoHashGrid_Imp.hpp
../src/geometry/GeoOctreeNode.hpp
../src/geometry/GeoOctreeNode_Imp.hpp
../src/geometry/GeoOctree.hpp
../src/geometry/GeoOctree_Imp.hpp
../src/geometry/GeoRtree.hpp
../src/geometry/GeoRtree_Imp.hpp
../src/geometry/GeoAdtree.hpp
../src/geometry/GeoAdtree_Imp.hpp
../src/geometry/GeoConvexHull.hpp
../src/geometry/GeoConvexHull_Imp.hpp
)

set(MESH_HEADERS
../src/mesh/MshNode.hpp
../src/mesh/MshFace.hpp
../src/mesh/MshFace_Imp.hpp
../src/mesh/MshElement.hpp
../src/mesh/MshElement_Imp.hpp
../src/mesh/MshMesh.hpp
../src/mesh/MshMesh_Imp.hpp
../src/mesh/MshCheckMesh.hpp
../src/mesh/MshCheckMesh_Imp.hpp
../src/mesh/MshMeshQuery.hpp
../src/mesh/MshMeshQuery_Imp.hpp
../src/mesh/MshBasicMesher.hpp
../src/mesh/MshBasicMesher_Imp.hpp
../src/mesh/MshExtrudedMesher.hpp
../src/mesh/MshExtrudedMesher_Imp.hpp
../src/mesh/MshTriangle.hpp
../src/mesh/MshTriangle_Imp.hpp
../src/mesh/MshQuadrilateral.hpp
../src/mesh/MshQuadrilateral_Imp.hpp
../src/mesh/MshTetrahedron.hpp
../src/mesh/MshTetrahedron_Imp.hpp
../src/mesh/MshTriangleMesher.hpp
../src/mesh/MshTriangleMesher_Imp.hpp
../src/mesh/MshQuadrilateralMesher.hpp
../src/mesh/MshQuadrilateralMesher_Imp.hpp
)

set(INT_HEADERS
../src/integration/IntGaussIntegration.hpp
../src/integration/IntHexahedron.hpp
../src/integration/IntHexahedron_Imp.hpp
../src/integration/IntQuadrilateral.hpp
../src/integration/IntQuadrilateral_Imp.hpp
../src/integration/IntTetrahedron.hpp
../src/integration/IntTetrahedron_Imp.hpp
../src/integration/IntTriangle.hpp
../src/integration/IntTriangle_Imp.hpp
../src/integration/IntTriangularPrism.hpp
../src/integration/IntTriangularPrism_Imp.hpp
)

set(BEM_HEADERS
../src/bem/BemTriangle.hpp
../src/bem/BemTriangle_Imp.hpp
)

set(FDM_HEADERS
../src/fdm/FdmGrid.hpp
../src/fdm/FdmGrid_Imp.hpp
)

set(FEM_HEADERS
../src/fem/FemElement.hpp
../src/fem/FemElement_Imp.hpp
../src/fem/FemEdge.hpp
../src/fem/FemEdge_Imp.hpp
../src/fem/FemFace.hpp
../src/fem/FemFace_Imp.hpp
../src/fem/FemSolid.hpp
../src/fem/FemBeam.hpp
../src/fem/FemBeam_Imp.hpp
../src/fem/FemTriangle.hpp
../src/fem/FemTriangle_Imp.hpp
../src/fem/FemQuadrilateral.hpp
../src/fem/FemQuadrilateral_Imp.hpp
../src/fem/FemTetrahedron.hpp
../src/fem/FemTetrahedron_Imp.hpp
../src/fem/FemHexahedron.hpp
../src/fem/FemHexahedron_Imp.hpp
../src/fem/FemTriangularPrism.hpp
../src/fem/FemTriangularPrism_Imp.hpp
)

set(FEM_THERMAL_HEADERS
../src/fem/thermal/FemThermalElement.hpp
../src/fem/thermal/FemThermalElement_Imp.hpp
../src/fem/thermal/FemLinearTemperatureBeam.hpp
../src/fem/thermal/FemLinearTemperatureBeam_Imp.hpp
../src/fem/thermal/FemLinearTemperatureTriangle.hpp
../src/fem/thermal/FemLinearTemperatureTriangle_Imp.hpp
../src/fem/thermal/FemLinearTemperatureQuadrilateral.hpp
../src/fem/thermal/FemLinearTemperatureQuadrilateral_Imp.hpp
../src/fem/thermal/FemLinearTemperatureTetrahedron.hpp
../src/fem/thermal/FemLinearTemperatureTetrahedron_Imp.hpp
../src/fem/thermal/FemLinearTemperatureTriangularPrism.hpp
../src/fem/thermal/FemLinearTemperatureTriangularPrism_Imp.hpp
../src/fem/thermal/FemLinearTemperatureHexahedron.hpp
../src/fem/thermal/FemLinearTemperatureHexahedron_Imp.hpp
)

set(FEM_STRUCTURAL_HEADERS
../src/fem/structural/FemStructuralElement.hpp
../src/fem/structural/FemStructuralElement_Imp.hpp
../src/fem/structural/FemConstantStrainTriangle.hpp
../src/fem/structural/FemConstantStrainTriangle_Imp.hpp
../src/fem/structural/FemConstantStrainTetrahedron.hpp
../src/fem/structural/FemConstantStrainTetrahedron_Imp.hpp
)

set(FEM_FLOW_HEADERS
../src/fem/flow/FemFlowElement.hpp
../src/fem/flow/FemFlowElement_Imp.hpp
../src/fem/flow/FemFlowTriangle.hpp
../src/fem/flow/FemFlowTriangle_Imp.hpp
../src/fem/flow/FemFlowTetrahedron.hpp
../src/fem/flow/FemFlowTetrahedron_Imp.hpp
)

set(FVM_HEADERS
../src/fvm/FvmNode.hpp
../src/fvm/FvmFace.hpp
../src/fvm/FvmFace_Imp.hpp
../src/fvm/FvmCell.hpp
../src/fvm/FvmCell_Imp.hpp
../src/fvm/FvmControlVolume.hpp
../src/fvm/FvmControlVolume_Imp.hpp
../src/fvm/FvmMesh.hpp
../src/fvm/FvmMesh_Imp.hpp
../src/fvm/FvmMeshSearch.hpp
../src/fvm/FvmMeshSearch_Imp.hpp
)

set(FVM_THERMAL_HEADERS
../src/fvm/thermal/FvmTemperatureSolver.hpp
../src/fvm/thermal/FvmTemperatureSolver_Imp.hpp
)

set(FVM_FLOW_HEADERS
../src/fvm/flow/FvmPisoSolver.hpp
../src/fvm/flow/FvmPisoSolver_Imp.hpp
../src/fvm/flow/FvmVofSolver.hpp
../src/fvm/flow/FvmVofSolver_Imp.hpp
)

set(SPH_HEADERS
../src/sph/SphKernel.hpp
../src/sph/SphKernel_Imp.hpp
../src/sph/SphCubicSpline.hpp
../src/sph/SphCubicSpline_Imp.hpp
../src/sph/SphGaussian.hpp
../src/sph/SphGaussian_Imp.hpp
../src/sph/SphQuintic.hpp
../src/sph/SphQuintic_Imp.hpp
../src/sph/SphSpiky.hpp
../src/sph/SphSpiky_Imp.hpp
../src/sph/SphConvex.hpp
../src/sph/SphConvex_Imp.hpp
../src/sph/SphParticles.hpp
../src/sph/SphParticles_Imp.hpp
)

set(SLE_HEADERS
../src/sle/SleSystem.hpp
../src/sle/SleSystem_Imp.hpp
)

set(PDE_HEADERS
../src/pde/PdeEquation.hpp
../src/pde/PdeEquation_Imp.hpp
../src/pde/PdeField.hpp
../src/pde/PdeField_Imp.hpp
../src/pde/PdeBoundaryCondition.hpp
../src/pde/PdeBoundaryCondition_Imp.hpp
)

set(PDE_FDM_HEADERS
../src/pde/fdm/FdmOperators.hpp
../src/pde/fdm/FdmOperators_Imp.hpp
)

set(PDE_FDM_GENERIC_HEADERS
../src/pde/fdm/generic/FdmGenericOperators_Imp.hpp
)

set(PDE_FDM_THERMAL_HEADERS
../src/pde/fdm/thermal/FdmThermalOperators_Imp.hpp
)

set(PDE_BEM_HEADERS
../src/pde/bem/BemOperators.hpp
../src/pde/bem/BemOperators_Imp.hpp
)

set(PDE_BEM_THERMAL_HEADERS
../src/pde/bem/thermal/BemThermalOperators_Imp.hpp
)

set(PDE_FEM_HEADERS
../src/pde/fem/FemOperators.hpp
../src/pde/fem/FemOperators_Imp.hpp
)

set(PDE_FEM_GENERIC_HEADERS
../src/pde/fem/generic/FemGenericOperators_Imp.hpp
)

set(PDE_FEM_THERMAL_HEADERS
../src/pde/fem/thermal/FemThermalOperators_Imp.hpp
)

set(PDE_FEM_STRUCTURAL_HEADERS
../src/pde/fem/structural/FemStructuralOperators_Imp.hpp
)

set(PDE_FEM_FLOW_HEADERS
../src/pde/fem/flow/FemFlowOperators_Imp.hpp
)

set(PDE_FVM_HEADERS
../src/pde/fvm/FvmOperators.hpp
../src/pde/fvm/FvmOperators_Imp.hpp
)

set(PDE_FVM_GENERIC_HEADERS
../src/pde/fvm/generic/FvmGenericOperators_Imp.hpp
)

set(PDE_FVM_THERMAL_HEADERS
../src/pde/fvm/thermal/FvmThermalOperators_Imp.hpp
)

set(PDE_FVM_FLOW_HEADERS
../src/pde/fvm/flow/FvmFlowOperators_Imp.hpp
)

set(MAT_HEADERS
../src/material/MatMaterial.hpp
../src/material/MatMaterial_Imp.hpp
)

set(ANA_HEADERS
../src/analytical/AnaFunction.hpp
../src/analytical/AnaFunction_Imp.hpp
../src/analytical/AnaTemperature.hpp
../src/analytical/AnaTemperature_Imp.hpp
)

set(POS_HEADERS
../src/../src/post/PosGmsh.hpp
../src/../src/post/PosGmsh_Imp.hpp
../src/post/PosVtk.hpp
../src/post/PosVtk_Imp.hpp
../src/post/PosGnuplot.hpp
../src/post/PosGnuplot_Imp.hpp
../src/post/PosQuickMesh.hpp
../src/post/PosQuickMesh_Imp.hpp
)

set(STL_HEADERS
../src/stl/StlUtils.hpp
../src/stl/StlUtils_Imp.hpp
../src/stl/StlFile.hpp
../src/stl/StlFile_Imp.hpp
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
source_group("Header Files\\pde\\fvm\\flow" FILES ${PDE_FVM_FLOW_HEADERS})
source_group("Header Files\\analytical" FILES ${ANA_HEADERS})
source_group("Header Files\\post" FILES ${POS_HEADERS})
source_group("Header Files\\stl" FILES ${STL_HEADERS})

source_group("Header Files" FILES ENigMA.hpp)
source_group("Source Files" FILES ENigMA.cpp)

include_directories(
    ../include
    ../src/common
    ../src/analytical
    ../src/bem
    ../src/geometry
    ../src/integration
    ../src/fdm
    ../src/fem
    ../src/fvm
    ../src/material
    ../src/mesh
    ../src/pde
    ../src/post
    ../src/sle
    ../src/stl
    ${Eigen3_DIR}
    ${ViennaCL_DIR}
    ${RTree_DIR}
    ${Exprtk_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(ENigMALib 
ENigMA.hpp
ENigMA.cpp
${GEOMETRY_HEADERS} 
${MESH_HEADERS} 
${MAT_HEADERS} 
${INT_HEADERS} 
${BEM_HEADERS} 
${FDM_HEADERS}
${FEM_HEADERS} 
${FEM_GENERIC_HEADERS} 
${FEM_THERMAL_HEADERS} 
${FEM_STRUCTURAL_HEADERS} 
${FEM_FLOW_HEADERS}
${FVM_HEADERS} 
${FVM_GENERIC_HEADERS} 
${FVM_THERMAL_HEADERS} 
${SPH_HEADERS}
${SLE_HEADERS} 
${PDE_HEADERS} 
${PDE_BEM_HEADERS}
${PDE_BEM_THERMAL_HEADERS}
${PDE_FDM_HEADERS}
${PDE_FDM_GENERIC_HEADERS}
${PDE_FDM_THERMAL_HEADERS}
${PDE_FEM_HEADERS} 
${PDE_FEM_GENERIC_HEADERS}
${PDE_FEM_THERMAL_HEADERS}
${PDE_FEM_STRUCTURAL_HEADERS}
${PDE_FEM_FLOW_HEADERS}
${PDE_FVM_HEADERS} 
${PDE_FVM_GENERIC_HEADERS}
${PDE_FVM_THERMAL_HEADERS}
${PDE_FVM_FLOW_HEADERS}
${ANA_HEADERS} 
${POS_HEADERS} 
${STL_HEADERS}
)


