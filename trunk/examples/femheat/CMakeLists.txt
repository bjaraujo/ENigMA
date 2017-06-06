
project(FemHeat)

set(FEMHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${FEMHEAT_SOURCES})

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
../../src/fem
../../src/bem
../../src/sle
../../src/pde
../../src/post
../../src/analytical
)

add_executable(FemHeat ${FEMHEAT_SOURCES})




