
project(FvmFlow)

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
    ../../src/fvm
    ../../src/sle
    ../../src/pde
    ../../src/post
    ../../src/analytical
)

# Set your files and resources here  
set(FVMFLOW_HEADERS FvmFlow.h)
set(FVMFLOW_SOURCES main.cpp FvmFlow.cpp)
set(FVMFLOW_UI FvmFlow.ui)
set(FVMFLOW_RESOURCES Icons/icons.qrc)

add_executable(FvmFlow ${FVMFLOW_HEADERS} ${FVMFLOW_SOURCES} ${FVMFLOW_UI} ${FVMFLOW_RESOURCES})

target_link_libraries(FvmFlow
  Qt5::Core
  Qt5::Widgets
  Qt5::Gui
  ${VTK_LIBRARIES}
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  add_definitions(/bigobj)
endif()



