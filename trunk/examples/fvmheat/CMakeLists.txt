
project(FvmHeat)

set(FVMHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${FVMHEAT_SOURCES})

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
  ../../src/fem
  ../../src/fdm
  ../../src/fvm
  ../../src/sle
  ../../src/pde
  ../../src/post
  ../../src/analytical
)

add_executable(FvmHeat ${FVMHEAT_SOURCES})

add_dependencies(FvmHeat eigen)


