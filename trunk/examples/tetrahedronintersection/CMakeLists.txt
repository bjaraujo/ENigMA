project(TetrahedronIntersection)

if(NOT VTK_BINARY_DIR)
  find_package(VTK)
  if(NOT VTK_DIR)
    message(FATAL_ERROR "Please set VTK_DIR.")
  endif(NOT VTK_DIR)
  include(${VTK_USE_FILE})
endif(NOT VTK_BINARY_DIR)

if(NOT Qt5_FOUND)
  find_package(Qt5Core)
  find_package(Qt5Gui)
  find_package(Qt5Widgets)
endif(NOT Qt5_FOUND)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Use the include path and library for Qt that is used by VTK.
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Qt5Core_INCLUDE_DIRS}
  ${Qt5Gui_INCLUDE_DIRS}
  ${Qt5Widgets_INCLUDE_DIRS}  
  ${Eigen3_DIR}
  ${RTree_DIR}
  ${Exprtk_DIR}
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

# Set your files and resources here  
set(TETRAHEDRONINTERSECTION_HEADERS TetrahedronIntersection.h)
set(TETRAHEDRONINTERSECTION_SOURCES main.cpp TetrahedronIntersection.cpp)
set(TETRAHEDRONINTERSECTION_UI TetrahedronIntersection.ui)
set(TETRAHEDRONINTERSECTION_RESOURCES Icons/icons.qrc)

add_executable(TetrahedronIntersection ${TETRAHEDRONINTERSECTION_HEADERS} ${TETRAHEDRONINTERSECTION_SOURCES} ${TETRAHEDRONINTERSECTION_UI} ${TETRAHEDRONINTERSECTION_RESOURCES})

target_link_libraries(TetrahedronIntersection
  Qt5::Core
  Qt5::Widgets
  Qt5::Gui
  ${VTK_LIBRARIES}
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()
