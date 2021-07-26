project(Mesh3D)

if(NOT VTK_BINARY_DIR)
  find_package(VTK)
  if(NOT VTK_DIR)
    message(FATAL_ERROR "Please set VTK_DIR.")
  endif(NOT VTK_DIR)
  #include(${VTK_USE_FILE})
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

# Set your files and resources here  
set(MESH3D_HEADERS Mesh3d.h)
set(MESH3D_SOURCES main.cpp Mesh3d.cpp)
set(MESH3D_UI Mesh3d.ui)
set(MESH3D_RESOURCES Icons/icons.qrc)

add_executable(Mesh3D ${MESH3D_HEADERS} ${MESH3D_SOURCES} ${MESH3D_UI} ${MESH3D_RESOURCES})

add_dependencies(Mesh3D eigen)

target_link_libraries(Mesh3D
  Qt5::Core
  Qt5::Widgets
  Qt5::Gui
  ${VTK_LIBRARIES}
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  add_definitions(/bigobj)
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()
