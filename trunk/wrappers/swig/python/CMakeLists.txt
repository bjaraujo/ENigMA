
project(ENigMA)

find_package(PythonLibs 3.0 REQUIRED)
find_package(PythonInterp 3.0 REQUIRED)

include_directories(
    ../../../src/common
    ../../../src/analytical
    ../../../src/bem
    ../../../src/csg
    ../../../src/geometry
    ../../../src/integration
    ../../../src/fdm
    ../../../src/fem
    ../../../src/fem/flow
    ../../../src/fem/structural
    ../../../src/fem/thermal
    ../../../src/fvm
    ../../../src/fvm/flow
    ../../../src/fvm/thermal
    ../../../src/sph
    ../../../src/lbm
    ../../../src/material
    ../../../src/mesh
    ../../../src/pde
    ../../../src/post
    ../../../src/sle
    ../../../src/stl
    ${Eigen3_DIR}
    ${Eigen3_DIR}/Eigen
    ${ViennaCL_DIR}
    ${RTree_DIR}
    ${Exprtk_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PYTHON_INCLUDE_PATH}
)

set(ENGIMA_INTERFACE
    ../ENigMA.i
)

set(CMAKE_SWIG_FLAGS "")

set_source_files_properties(${ENGIMA_INTERFACE} PROPERTIES CPLUSPLUS ON)
set_source_files_properties(${ENGIMA_INTERFACE} PROPERTIES SWIG_FLAGS "-debug-tmused")

swig_add_module(ENigMA python ${ENGIMA_INTERFACE})
swig_link_libraries(ENigMA ${PYTHON_LIBRARIES})

if(WIN32)
add_definitions(/bigobj)
endif(WIN32)

add_custom_command(TARGET _ENigMA POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ENigMA.py $<TARGET_FILE_DIR:_ENigMA>/ENigMA.py)

