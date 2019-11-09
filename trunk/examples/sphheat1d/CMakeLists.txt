
project(SphHeat1D)

if(NOT Qt5_FOUND)
  find_package(Qt5Core)
  find_package(Qt5Gui)
  find_package(Qt5Widgets)
  find_package(Qt5PrintSupport)
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
  ${VIENNACL_DIR}
  ${RTREE_DIR}
  ${EXPRTK_DIR}
  ../../src/common
  ../../src/geometry
  ../../src/integration
  ../../src/mesh
  ../../src/material
  ../../src/bem
  ../../src/fdm
  ../../src/fem
  ../../src/fvm
  ../../src/sle
  ../../src/lbm
  ../../src/sph
  ../../src/pde
  ../../src/post
  ../../src/analytical
)

set(SPH_HEAT1D_HEADERS MainWindow.h QCustomPlot.h)
set(SPH_HEAT1D_SOURCES main.cpp MainWindow.cpp QCustomPlot.cpp)
set(SPH_HEAT1D_UI MainWindow.ui)

add_executable(SphHeat1D ${SPH_HEAT1D_HEADERS} ${SPH_HEAT1D_SOURCES} ${SPH_HEAT1D_UI})

target_link_libraries(SphHeat1D
  Qt5::Core
  Qt5::Widgets
  Qt5::Gui
  Qt5::PrintSupport
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  add_definitions(/bigobj)
endif()



