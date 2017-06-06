
project(FemFlow)

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
  ../../src/material
  ../../src/bem
  ../../src/fem
  ../../src/fem/flow
  ../../src/fvm
  ../../src/sle
  ../../src/pde
  ../../src/post
  ../../src/analytical
)

# Set your files and resources here  
set(FEMFLOW_HEADERS FemFlow.h)
set(FEMFLOW_SOURCES main.cpp FemFlow.cpp)
set(FEMFLOW_UI FemFlow.ui)
set(FEMFLOW_RESOURCES Icons/icons.qrc)

add_executable(FemFlow ${FEMFLOW_HEADERS} ${FEMFLOW_SOURCES} ${FEMFLOW_UI} ${FEMFLOW_RESOURCES})

target_link_libraries(FemFlow
  Qt5::Core
  Qt5::Widgets
  Qt5::Gui
  ${VTK_LIBRARIES}
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  add_definitions(/bigobj)
endif()



