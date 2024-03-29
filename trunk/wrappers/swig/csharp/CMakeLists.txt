
project(ENigMAcs)

find_package(CSharp REQUIRED)
include(${CSHARP_USE_FILE})

include_directories(
  ${EIGEN_DIR}
  ${EIGEN_DIR}/Eigen
  ${VIENNACL_DIR}
  ${RTREE_DIR}
  ${EXPRTK_DIR}
  ../../../src/common
  ../../../src/analytical
  ../../../src/bem
  ../../../src/csg
  ../../../src/geometry
  ../../../src/integration
  ../../../src/fdm
  ../../../src/fem
  ../../../src/fvm
  ../../../src/fvm/flow
  ../../../src/fvm/thermal
  ../../../src/sph
  ../../../src/lbm
  ../../../src/material
  ../../../src/mesh
  ../../../src/pde
  ../../../src/solvers
  ../../../src/post
  ../../../src/sle
  ../../../src/stl
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PYTHON_INCLUDE_PATH}
)

set(ENGIMA_INTERFACE
  ../ENigMA.i
)

#set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(${ENGIMA_INTERFACE} PROPERTIES CPLUSPLUS ON)

swig_add_library(ENigMAcs LANGUAGE csharp SOURCES ${ENGIMA_INTERFACE})

add_dependencies(ENigMAcs eigen)

if(WIN32)
  add_definitions(/bigobj)
endif()

install(TARGETS ENigMAcs DESTINATION ENigMAcs)
install(FILES "${CMAKE_SOURCE_DIR}/../LICENSE.txt" DESTINATION ENigMAcs)
install(FILES "${CMAKE_SOURCE_DIR}/../README.md" DESTINATION ENigMAcs)
install(FILES "${CMAKE_SOURCE_DIR}/src/version.h" DESTINATION ENigMAcs)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" DESTINATION ENigMAcs)
