project(ClipPolyhedron)

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
  ${EIGEN_DIR}
  ${RTREE_DIR}
  ${EXPRTK_DIR}
  ../../src/common
  ../../src/geometry
)

# Set your files and resources here  
set(CLIPPOLYHEDRON_HEADERS ClipPolyhedron.h)
set(CLIPPOLYHEDRON_SOURCES main.cpp ClipPolyhedron.cpp)
set(CLIPPOLYHEDRON_UI ClipPolyhedron.ui)
set(CLIPPOLYHEDRON_RESOURCES Icons/icons.qrc)

add_executable(ClipPolyhedron ${CLIPPOLYHEDRON_HEADERS} ${CLIPPOLYHEDRON_SOURCES} ${CLIPPOLYHEDRON_UI})

target_link_libraries(ClipPolyhedron
  Qt5::Core
  Qt5::Widgets
  Qt5::Gui
  ${VTK_LIBRARIES}
)


