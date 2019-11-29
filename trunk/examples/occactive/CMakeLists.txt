
project(OccActive)

find_package(OpenCASCADE REQUIRED COMPONENTS TKernel TKMath TKService TKV3d TKOpenGl TKBRep TKIGES TKSTL TKVRML TKSTEP TKSTEPAttr TKSTEP209 TKSTEPBase TKGeomBase TKGeomAlgo TKG3d TKG2d TKXSBase TKShHealing TKHLR TKTopAlgo TKMesh TKPrim TKCDF TKBool TKBO TKFillet TKOffset)

if(NOT Qt5_FOUND)
  find_package(Qt5Core)
  find_package(Qt5Gui)
  find_package(Qt5Widgets)
  find_package(Qt5PrintSupport)
endif(NOT Qt5_FOUND)

include(${QT_USE_FILE})

include_directories(
    ../../src/common
    ../../src/analytical
    ../../src/geometry
    ../../src/mesh
    ../../src/material
    ../../src/integration
    ../../src/fvm
    ../../src/sle
    ../../src/bem
    ../../src/pde
    ../../src/post
    ../../src/stl
    ${EIGEN_DIR}
    ${VIENNACL_DIR}
    ${RTREE_DIR}
    ${EXPRTK_DIR}
    ${OpenCASCADE_INCLUDE_DIR}
    ${QT_INCLUDE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set your files and resources here  
set(OCC_ACTIVEQT_HEADERS src/qt/Manager.h)
set(OCC_ACTIVEQT_SOURCES )

set(OCC_ACTIVEQT_MOC_HEADERS src/qt/MainWindow.h)
set(OCC_ACTIVEQT_MOC_SOURCES src/qt/main.cpp src/qt/MainWindow.cpp)
set(OCC_ACTIVEQT_UI src/qt/MainWindow.ui)
set(OCC_ACTIVEQT_RESOURCES src/qt/MainWindow.qrc)

set(OCC_ACTIVEQT_MOC_OCC_HEADERS src/qt/occ/OccWidget.h src/qt/occ/OccWindow.h src/qt/occ/OccManager.h src/qt/occ/OccStructureView.h src/qt/occ/OccStructureModel.h src/qt/occ/OccPropertyView.h src/qt/occ/OccPropertyModel.h)
set(OCC_ACTIVEQT_MOC_OCC_SOURCES src/qt/occ/OccWidget.cpp src/qt/occ/OccWindow.cpp src/qt/occ/OccManager.cpp src/qt/occ/OccStructureView.cpp src/qt/occ/OccStructureModel.cpp src/qt/occ/OccPropertyView.cpp src/qt/occ/OccPropertyModel.cpp)

set(OCC_ACTIVEOCC_HEADERS src/occ/OccTranslator.h src/occ/OccMesher.h)
set(OCC_ACTIVEOCC_SOURCES src/occ/OccTranslator.cpp src/occ/OccMesher.cpp)

# The rest should just work (sure...)
QT4_WRAP_UI(UISrcs ${OCC_ACTIVEQT_UI})
QT4_WRAP_CPP(MOCSrcs ${OCC_ACTIVEQT_MOC_HEADERS} ${OCC_ACTIVEQT_MOC_OCC_HEADERS})
QT4_ADD_RESOURCES(ResourceSrcs ${OCC_ACTIVEQT_RESOURCES})

source_group("Header Files\\qt" FILES
    ${OCC_ACTIVEQT_HEADERS}
    ${OCC_ACTIVEQT_MOC_HEADERS}
)

source_group("Header Files\\occ" FILES
    ${OCC_ACTIVEOCC_HEADERS}
)

source_group("Header Files\\qt\\occ" FILES
    ${OCC_ACTIVEQT_MOC_OCC_HEADERS}
)

source_group("Source Files\\qt" FILES
    ${OCC_ACTIVEQT_SOURCES}
    ${OCC_ACTIVEQT_MOC_SOURCES}
)

source_group("Source Files\\occ" FILES
    ${OCC_ACTIVEOCC_SOURCES}
)

source_group("Source Files\\qt\\occ" FILES
    ${OCC_ACTIVEQT_MOC_OCC_SOURCES}
)

source_group("Resource Files\\qt" FILES
    ${OCC_ACTIVEQT_UI}
    ${OCC_ACTIVEQT_RESOURCES}
    ${EXE_ICON}
)

source_group("Generated Files\\qt" FILES
    ${UISrcs}
    ${MOCSrcs}
    ${ResourceSrcs}
    ${RCS_SOURCES}
)

add_definitions(-DQT_GUI_LIBS -DQT_CORE_LIB -DNOMINMAX)
add_definitions(-DCSFDB)

if(WIN32)
add_definitions(-DWNT)
add_definitions(/bigobj)
endif(WIN32)

set_source_files_properties(${OCC_ACTIVEQT_SOURCES} ${OCC_ACTIVEOCC_SOURCES} ${OCC_ACTIVEQT_OCC_SOURCES} PROPERTIES OBJECT_DEPENDS "${UISrcs}")

set(OCC_ACTIVEHEADERS ${OCC_ACTIVEQT_HEADERS} ${OCC_ACTIVEQT_MOC_HEADERS} ${OCC_ACTIVEOCC_HEADERS} ${OCC_ACTIVEQT_MOC_OCC_HEADERS} ${QT_SOLUTION_HEADERS})
set(OCC_ACTIVESOURCES ${OCC_ACTIVEQT_SOURCES} ${OCC_ACTIVEQT_MOC_SOURCES} ${OCC_ACTIVEOCC_SOURCES} ${OCC_ACTIVEQT_MOC_OCC_SOURCES} ${QT_SOLUTION_SOURCES})

IF (USE_NETGEN)
include_directories(${NETGEN_DIR}/nglib)
add_definitions(-DUSE_NETGEN)
link_directories(${NETGEN_LIBDIR})
endif(USE_NETGEN)

add_executable(OccActive MACOSX_BUNDLE ${OCC_ACTIVEHEADERS} ${OCC_ACTIVESOURCES} ${UISrcs} ${MOCSrcs} ${ResourceSrcs})

add_dependencies(OccActive eigen)

target_link_libraries(OccActive
    ${QT_LIBRARIES} 
    ${OpenCASCADE_LIBRARIES}
)

IF (USE_NETGEN)
target_link_libraries(OccActive nglib)
endif(USE_NETGEN)



