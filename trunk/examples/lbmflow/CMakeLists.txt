
project(LbmFlow)

set(LBMFLOW_SOURCES
main.cpp
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

source_group("Source Files" FILES ${LBMFLOW_SOURCES})

include_directories(
${Eigen3_DIR}
${RTree_DIR}
${Exprtk_DIR}
../../src/common
../../src/geometry
../../src/integration
../../src/mesh
../../src/material
../../src/sle
../../src/lbm
../../src/post
../../src/solvers
../../src/analytical
${OPENGL_INCLUDE_DIRS}
${GLUT_INCLUDE_DIRS}
)

add_executable(LbmFlow ${LBMFLOW_SOURCES})

target_link_libraries(LbmFlow
${OPENGL_LIBRARIES}
${GLUT_LIBRARY}
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_definitions(/bigobj)
endif()



