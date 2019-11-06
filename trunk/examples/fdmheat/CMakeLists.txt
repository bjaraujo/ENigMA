
project(FdmHeat)

set(FDMHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${FDMHEAT_SOURCES})

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
  ../../src/fvm
  ../../src/fdm
  ../../src/sle
  ../../src/pde
  ../../src/post
  ../../src/analytical
)

add_executable(FdmHeat ${FDMHEAT_SOURCES})




