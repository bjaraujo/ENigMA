
project(UnitTests)

set(BUILD_SHARED_LIBS ON)

set(TEST_GEOMETRY_SOURCES
TestGeoNormal.cpp
TestGeoCoordinate.cpp
TestGeoBoundingBox.cpp
TestGeoTransform2D.cpp
TestGeoVector.cpp
TestGeoPlane.cpp
TestGeoLine.cpp
TestGeoLineList.cpp
TestGeoPolygon.cpp
TestGeoPolyline.cpp
TestGeoTriangle.cpp
TestGeoQuadrilateral.cpp
TestGeoTetrahedron.cpp
TestGeoHexahedron.cpp
TestGeoTriangularPrism.cpp
TestGeoPolyhedron.cpp
TestGeoCircle.cpp
TestGeoSphere.cpp
TestGeoHashGrid.cpp
TestGeoOctree.cpp
TestGeoRtree.cpp
TestGeoAdtree.cpp
TestGeoConvexHull.cpp
)

set(TEST_MESH_SOURCES
TestMshNode.cpp
TestMshElement.cpp
TestMshMesh.cpp
TestMshMeshQuery.cpp
TestMshBasicMesher.cpp
TestMshExtrudedMesher.cpp
TestMshTriangleMesher.cpp
TestMshQuadrilateralMesher.cpp
TestMshTetrahedronMesher.cpp
)

set(TEST_BEM_SOURCES
TestBemTriangle.cpp
)

set(TEST_FEM_SOURCES
TestFemBeam.cpp
TestFemTriangle.cpp
TestFemQuadrilateral.cpp
TestFemTetrahedron.cpp
TestFemTriangularPrism.cpp
TestFemHexahedron.cpp
TestFemCbs.cpp
)

set(TEST_FVM_SOURCES
TestFvmFace.cpp
TestFvmControlVolume.cpp
TestFvmMesh.cpp
TestFvmPiso.cpp
)

set(TEST_SPH_SOURCES
TestSphCubicSpline.cpp
TestSphGaussian.cpp
TestSphQuintic.cpp
TestSphSpiky.cpp
TestSphConvex.cpp
)

set(TEST_PDE_SOURCES
TestPdeEquation.cpp
TestPdeField.cpp
)

set(TEST_ANA_SOURCES
TestAnaFunction.cpp
TestAnaTemperature.cpp
)

set(TEST_POS_SOURCES
TestPosGmsh.cpp
TestPosVtk.cpp
TestPosGnuplot.cpp
)

set(TEST_STL_SOURCES
TestStlUtils.cpp
)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests/googletest-download"
)
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests/googletest-download"
)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds the following targets:
# gtest, gtest_main, gmock and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build"
)

enable_testing()

source_group("Source Files\\geo" FILES ${TEST_GEOMETRY_SOURCES})
source_group("Source Files\\mesh" FILES ${TEST_MESH_SOURCES})
source_group("Source Files\\bem" FILES ${TEST_BEM_SOURCES})
source_group("Source Files\\fem" FILES ${TEST_FEM_SOURCES})
source_group("Source Files\\fvm" FILES ${TEST_FVM_SOURCES})
source_group("Source Files\\sph" FILES ${TEST_SPH_SOURCES})
source_group("Source Files\\pde" FILES ${TEST_PDE_SOURCES})
source_group("Source Files\\ana" FILES ${TEST_ANA_SOURCES})
source_group("Source Files\\post" FILES ${TEST_POS_SOURCES})
source_group("Source Files\\stl" FILES ${TEST_STL_SOURCES})

include_directories(
  ${EIGEN_DIR}
  ${VIENNACL_DIR}
  ${RTREE_DIR}
  ${EXPRTK_DIR}
  ../src/common
  ../src/geometry
  ../src/mesh
  ../src/material
  ../src/integration
  ../src/bem
  ../src/fdm
  ../src/fem
  ../src/fem/flow
  ../src/fvm
  ../src/sph
  ../src/sle
  ../src/pde
  ../src/solvers
  ../src/analytical
  ../src/post
  ../src/stl
)

add_executable(UnitTests ${TEST_GEOMETRY_SOURCES} ${TEST_MESH_SOURCES} ${TEST_BEM_SOURCES} ${TEST_FEM_SOURCES} ${TEST_CVFEM_SOURCES} ${TEST_FVM_SOURCES} ${TEST_SPH_SOURCES} ${TEST_PDE_SOURCES} ${TEST_ANA_SOURCES} ${TEST_POS_SOURCES} ${TEST_STL_SOURCES})

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_definitions(/bigobj)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

if(WIN32)
    target_link_libraries(UnitTests
    gtest
    gtest_main
)
else(WIN32)
    target_link_libraries(UnitTests
    gtest
    gtest_main
    pthread
)
endif(WIN32)

