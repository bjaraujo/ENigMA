
project(SphHeat)

set(SPHHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${SPHHEAT_SOURCES})

include_directories(
${Eigen3_DIR}
${ViennaCL_DIR}
${RTree_DIR}
${Exprtk_DIR}
  ../../src/common
  ../../src/geometry
  ../../src/integration
  ../../src/mesh
  ../../src/material
  ../../src/sph
  ../../src/bem
  ../../src/fdm
  ../../src/fem
  ../../src/fvm
  ../../src/sle
  ../../src/pde
  ../../src/post
  ../../src/analytical
)

add_executable(SphHeat ${SPHHEAT_SOURCES})




