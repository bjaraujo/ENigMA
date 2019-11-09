
project(LbmFlow)

set(LBMFLOW_SOURCES
main.cpp
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

source_group("Source Files" FILES ${LBMFLOW_SOURCES})

include_directories(
${EIGEN_DIR}
${RTREE_DIR}
${EXPRTK_DIR}
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
${GLUT_INCLUDE_DIR}
)

add_executable(LbmFlow ${LBMFLOW_SOURCES})

target_link_libraries(LbmFlow
${OPENGL_LIBRARIES}
${GLUT_LIBRARIES}
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_definitions(/bigobj)
endif()



