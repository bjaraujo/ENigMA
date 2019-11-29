project(Mesher)

set(MESHER_SOURCES
main.cpp
)

source_group("Source Files" FILES ${MESHER_SOURCES})

include_directories(
${EIGEN_DIR}
${RTREE_DIR}
${EXPRTK_DIR}
../../src/common
../../src/geometry
../../src/csg
../../src/integration
../../src/mesh
../../src/material
../../src/fvm
../../src/sle
../../src/pde
../../src/post
../../src/analytical
../../src/stl
)

add_executable(Mesher ${MESHER_SOURCES})

add_dependencies(Mesher eigen)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_definitions(/bigobj)
endif()

