
#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.8)

if(POLICY CMP0020)
  # Accept new policy related to automatically linking Qt apps to qtmain.lib on Windows
  cmake_policy(SET CMP0020 NEW)
endif()

project(ENigMA CXX)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)

#Configure options
set(EXPRTK_DIR "${CMAKE_SOURCE_DIR}/dependencies/exprtk/exprtk-2019")
set(RTREE_DIR "${CMAKE_SOURCE_DIR}/dependencies/RTree/RTree-2011")

set(QT_MOC_EXECUTABLE "$ENV{QTDIR}/bin/moc.exe" CACHE FILEPATH "")
set(QT_UIC_EXECUTABLE "$ENV{QTDIR}/bin/uic.exe" CACHE FILEPATH "")
set(QT_QMAKE_EXECUTABLE "$ENV{QTDIR}/bin/qmake.exe" CACHE FILEPATH "")
set(VTK_DIR $ENV{VTK_DIR} CACHE FILEPATH "")
set(CASROOT $ENV{CASROOT} CACHE FILEPATH "")

set(USE_VIENNACL OFF CACHE BOOL "Use ViennaCL")

set(ENIGMA_BUILD_CONSOLE_EXAMPLES OFF CACHE BOOL "Build console examples")
set(ENIGMA_BUILD_OPENGL_EXAMPLES OFF CACHE BOOL "Build OpenGL examples")
set(ENIGMA_BUILD_QT_EXAMPLES OFF CACHE BOOL "Build QT examples")
set(ENIGMA_BUILD_VTK_EXAMPLES OFF CACHE BOOL "Build VTK examples")
set(ENIGMA_BUILD_OCC_EXAMPLES OFF CACHE BOOL "Build OpenCASCADE examples")

set(ENIGMA_BUILD_UNIT_TESTS ON CACHE BOOL "Build unit tests")

set(ENIGMA_BUILD_WRAPPERS_SWIG OFF CACHE BOOL "Build SWIG wrapper")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY AUTOGEN_TARGETS_FOLDER automoc)

if(MSVC)
    add_definitions("/MP")
endif()

if(USE_VIENNACL)
    add_definitions(-DUSE_VIENNACL)
endif()

include(GetEigen)

if(ENIGMA_BUILD_UNIT_TESTS)
    add_subdirectory(tests)

    if (CMAKE_BUILD_TYPE STREQUAL "Coverage")
        set(Coverage_NAME UnitTests_coverage)
        include(CodeCoverage)
        setup_target_for_coverage_lcov(UnitTests_coverage UnitTests coverage)

        SET(CMAKE_CXX_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
        SET(CMAKE_C_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    endif()
endif()

add_subdirectory(wrappers)
add_subdirectory(examples)

