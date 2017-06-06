
project(FemIntro)

set(FEMINTRO_SOURCES
main.cpp
)

source_group("Source Files" FILES ${FEMINTRO_SOURCES})

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
../../src/sle
../../src/pde
../../src/analytical
../../src/material
../../src/post
)

add_executable(FemIntro ${FEMINTRO_SOURCES})



