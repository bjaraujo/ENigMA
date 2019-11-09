
project(SphFlow)


set(SPHFLOW_HEADERS
)

set(SPHFLOW_SOURCES
main.cpp
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

source_group("Header Files" FILES ${SPHFLOW_HEADERS})
source_group("Source Files" FILES ${SPHFLOW_SOURCES})

include_directories(
${EIGEN_DIR}
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
../../src/sle
../../src/pde
../../src/post
../../src/analytical
${OPENGL_INCLUDE_DIRS}
${GLUT_INCLUDE_DIR}
)

add_executable(SphFlow ${SPHFLOW_HEADERS} ${SPHFLOW_SOURCES})

target_link_libraries(SphFlow 
${OPENGL_LIBRARIES}
${GLUT_LIBRARIES}
)
