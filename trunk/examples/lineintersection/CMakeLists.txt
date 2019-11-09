
project(LineIntersection)

if(NOT Qt5_FOUND)
  find_package(Qt5Core)
  find_package(Qt5Gui)
  find_package(Qt5Widgets)
endif(NOT Qt5_FOUND)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

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
  ../../src/integration
  ../../src/mesh
  ../../src/fdm
  ../../src/fem
  ../../src/bem
  ../../src/fvm
  ../../src/sle
  ../../src/pde
  ../../src/analytical
  ../../src/material
  ../../src/post
)

set(LINEINTERSECTION_HEADERS MainWindow.h MyGraphicsView.h)
set(LINEINTERSECTION_SOURCES main.cpp MainWindow.cpp MyGraphicsView.cpp)
set(LINEINTERSECTION_UI MainWindow.ui)

add_executable(LineIntersection ${LINEINTERSECTION_HEADERS} ${LINEINTERSECTION_SOURCES} ${LINEINTERSECTION_UI})

target_link_libraries(LineIntersection
  Qt5::Core
  Qt5::Widgets
  Qt5::Gui
)


