
project(SphHeat)

set(SPHHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${SPHHEAT_SOURCES})

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
  ../../src/sph
  ../../src/bem
  ../../src/fdm
  ../../src/fem
  ../../src/fvm
  ../../src/sle
  ../../src/pde
  ../../src/post
  ../../src/analytical
  ../../src/solvers
)

add_executable(SphHeat ${SPHHEAT_SOURCES})

add_dependencies(SphHeat eigen)


