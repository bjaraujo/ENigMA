
project(FvmHeat)

set(FVMHEAT_SOURCES
main.cpp
)

source_group("Source Files" FILES ${FVMHEAT_SOURCES})

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
../../src/fvm
../../src/sle
../../src/pde
../../src/post
../../src/analytical
)

add_executable(FvmHeat ${FVMHEAT_SOURCES})



