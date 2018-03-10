
project(BemHeat)

set(BEMHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${BEMHEAT_SOURCES})

include_directories(
  ${Eigen3_DIR}
  ${ViennaCL_DIR}
  ${RTree_DIR}
  ${Exprtk_DIR}
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

add_executable(BemHeat ${BEMHEAT_SOURCES})

