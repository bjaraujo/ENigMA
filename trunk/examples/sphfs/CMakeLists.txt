
project(SphFS)


set(SPHFS_HEADERS
)

set(SPHFS_SOURCES
main.cpp
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

source_group("Header Files" FILES ${SPHFS_HEADERS})
source_group("Source Files" FILES ${SPHFS_SOURCES})

include_directories(
${Eigen3_DIR}
${RTree_DIR}
${Exprtk_DIR}
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
${GLUT_INCLUDE_DIRS}
)

add_executable(SphFS ${SPHFS_HEADERS} ${SPHFS_SOURCES})

target_link_libraries(SphFS 
${OPENGL_LIBRARIES}
${GLUT_LIBRARY}
)
