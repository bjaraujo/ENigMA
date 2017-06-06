
project(_ENigMAcs)

find_package(CSharp REQUIRED)
include(${CSHARP_USE_FILE})

include_directories(
    ../../../src/common
    ../../../src/analytical
    ../../../src/bem
    ../../../src/csg
    ../../../src/geometry
    ../../../src/integration
    ../../../src/fdm
    ../../../src/fem
    ../../../src/fvm
    ../../../src/sph
    ../../../src/lbm
    ../../../src/material
    ../../../src/mesh
    ../../../src/pde
    ../../../src/post
    ../../../src/sle
    ../../../src/stl
    ${Eigen3_DIR}
    ${Eigen3_DIR}/Eigen
    ${ViennaCL_DIR}
    ${RTree_DIR}
    ${Exprtk_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PYTHON_INCLUDE_PATH}
)

set(ENGIMA_INTERFACE
    ../ENigMA.i
)

set(CMAKE_SWIG_FLAGS "")

set_source_files_properties(${ENGIMA_INTERFACE} PROPERTIES CPLUSPLUS ON)

SWIG_ADD_MODULE(_ENigMAcs csharp ${ENGIMA_INTERFACE})

if(WIN32)
add_definitions(/bigobj)
endif(WIN32)
