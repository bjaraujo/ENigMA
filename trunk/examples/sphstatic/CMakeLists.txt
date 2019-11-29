
project(SphStatic)


set(SPHSTATIC_HEADERS
)

set(SPHSTATIC_SOURCES
main.cpp
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

source_group("Header Files" FILES ${SPHSTATIC_HEADERS})
source_group("Source Files" FILES ${SPHSTATIC_SOURCES})

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

add_executable(SphStatic ${SPHSTATIC_HEADERS} ${SPHSTATIC_SOURCES})

add_dependencies(SphStatic eigen)

target_link_libraries(SphStatic 
${OPENGL_LIBRARIES}
${GLUT_LIBRARIES}
)
