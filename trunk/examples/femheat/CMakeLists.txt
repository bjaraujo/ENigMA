
project(FemHeat)

set(FEMHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${FEMHEAT_SOURCES})

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
../../src/fdm
../../src/fem
../../src/fvm
../../src/bem
../../src/sle
../../src/pde
../../src/post
../../src/analytical
)

add_executable(FemHeat ${FEMHEAT_SOURCES})




