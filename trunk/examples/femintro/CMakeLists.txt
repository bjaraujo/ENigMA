
project(FemIntro)

set(FEMINTRO_SOURCES
main.cpp
)

source_group("Source Files" FILES ${FEMINTRO_SOURCES})

include_directories(
${EIGEN_DIR}
${VIENNACL_DIR}
${RTREE_DIR}
${EXPRTK_DIR}
  ../../src/common
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
  ../../src/analytical
  ../../src/material
  ../../src/post
)

add_executable(FemIntro ${FEMINTRO_SOURCES})



