
project(Finance)

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
  ${Eigen3_DIR}
  ${ViennaCL_DIR}
  ${RTree_DIR}
  ${Exprtk_DIR}
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

set(FINANCE_HEADERS MainWindow.h QCustomPlot.h BlackScholes.h)
set(FINANCE_SOURCES main.cpp MainWindow.cpp QCustomPlot.cpp BlackScholes.cpp)
set(FINANCE_UI MainWindow.ui)

add_executable(Finance ${FINANCE_HEADERS} ${FINANCE_SOURCES} ${FINANCE_UI})

target_link_libraries(Finance
  Qt5::Core
  Qt5::Widgets
  Qt5::Gui
  Qt5::PrintSupport
)


