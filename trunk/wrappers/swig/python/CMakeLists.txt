
project(ENigMA)

find_package(PythonLibs 3.0 REQUIRED)
find_package(PythonInterp 3.0 REQUIRED)

include_directories(
  ${EIGEN_DIR}
  ${EIGEN_DIR}/Eigen
  ${VIENNACL_DIR}
  ${RTREE_DIR}
  ${EXPRTK_DIR}
  ../../../src/common
  ../../../src/analytical
  ../../../src/bem
  ../../../src/csg
  ../../../src/geometry
  ../../../src/integration
  ../../../src/fdm
  ../../../src/fem
  ../../../src/fvm
  ../../../src/fvm/flow
  ../../../src/fvm/thermal
  ../../../src/sph
  ../../../src/lbm
  ../../../src/material
  ../../../src/mesh
  ../../../src/pde
  ../../../src/solvers
  ../../../src/post
  ../../../src/sle
  ../../../src/stl
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PYTHON_INCLUDE_PATH}
)

set(ENGIMA_INTERFACE
  ../ENigMA.i
)

#set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(${ENGIMA_INTERFACE} PROPERTIES CPLUSPLUS ON)

swig_add_library(ENigMA LANGUAGE python SOURCES ${ENGIMA_INTERFACE})
swig_link_libraries(ENigMA ${PYTHON_LIBRARIES})

add_dependencies(_ENigMA eigen)

if(WIN32)
  add_definitions(/bigobj)
endif(WIN32)

add_custom_command(TARGET _ENigMA POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ENigMA.py $<TARGET_FILE_DIR:_ENigMA>/ENigMA.py)

install(TARGETS _ENigMA DESTINATION ENigMApy/ENigMA)
install(FILES $<TARGET_FILE_DIR:_ENigMA>/ENigMA.py DESTINATION ENigMApy/ENigMA)
install(FILES "${CMAKE_SOURCE_DIR}/../LICENSE.txt" DESTINATION ENigMApy)
install(FILES "${CMAKE_SOURCE_DIR}/../README.md" DESTINATION ENigMApy)
install(FILES "${CMAKE_SOURCE_DIR}/src/version.h" DESTINATION ENigMApy)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in" DESTINATION ENigMApy)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" DESTINATION ENigMApy)



