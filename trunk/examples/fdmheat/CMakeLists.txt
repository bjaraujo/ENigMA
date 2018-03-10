
project(FdmHeat)

set(FDMHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${FDMHEAT_SOURCES})

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




