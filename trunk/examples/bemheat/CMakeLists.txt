
project(BemHeat)

set(BEMHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${BEMHEAT_SOURCES})

include_directories(
  ${EIGEN_DIR}
  ${VIENNACL_DIR}
  ${RTREE_DIR}
  ${EXPRTK_DIR}
  ../../src/common
  ../../src/analytical
  ../../src/geometry
  ../../src/integration
  ../../src/mesh
  ../../src/material
  ../../src/bem
  ../../src/fdm
  ../../src/fem
  ../../src/fvm
  ../../src/sle
  ../../src/pde
  ../../src/post
)

if(WIN32)
add_definitions(/bigobj)
endif(WIN32)

add_executable(BemHeat ${BEMHEAT_SOURCES})

add_dependencies(BemHeat eigen)
