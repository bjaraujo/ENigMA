project(OccActive)

include(${QT_USE_FILE})

include_directories(
    ${EIGEN_DIR}
    ${VIENNACL_DIR}
    ${RTREE_DIR}
    ${EXPRTK_DIR}
    ${OpenCASCADE_INCLUDE_DIR}
    ${QT_INCLUDE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
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
)

# Set your files and resources here  
set(OCC_ACTIVEQT_HEADERS qt/Manager.h)
set(OCC_ACTIVEQT_SOURCES )

set(OCC_ACTIVEQT_MOC_HEADERS qt/MainWindow.h)
set(OCC_ACTIVEQT_MOC_SOURCES qt/main.cpp qt/MainWindow.cpp)
set(OCC_ACTIVEQT_UI qt/MainWindow.ui)
set(OCC_ACTIVEQT_RESOURCES qt/MainWindow.qrc)

set(OCC_ACTIVEQT_MOC_OCC_HEADERS qt/occ/OccWidget.h qt/occ/OccWindow.h qt/occ/OccManager.h qt/occ/OccStructureView.h qt/occ/OccStructureModel.h qt/occ/OccPropertyView.h qt/occ/OccPropertyModel.h)
set(OCC_ACTIVEQT_MOC_OCC_SOURCES qt/occ/OccWidget.cpp qt/occ/OccWindow.cpp qt/occ/OccManager.cpp qt/occ/OccStructureView.cpp qt/occ/OccStructureModel.cpp qt/occ/OccPropertyView.cpp qt/occ/OccPropertyModel.cpp)

set(OCC_ACTIVEOCC_HEADERS occ/OccTranslator.h occ/OccMesher.h)
set(OCC_ACTIVEOCC_SOURCES occ/OccTranslator.cpp occ/OccMesher.cpp)

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

target_link_libraries(OccActive
    ${QT_LIBRARIES} 
    ${OpenCASCADE_LIBRARIES}
)

IF (USE_NETGEN)
target_link_libraries(OccActive nglib)
endif(USE_NETGEN)



