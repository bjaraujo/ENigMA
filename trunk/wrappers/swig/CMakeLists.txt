
project(SwigWrappers)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(CMAKE_SWIG_FLAGS "")

set(WRAP_SWIG_PYTHON OFF CACHE BOOL "Build SWIG python interface")

IF (WRAP_SWIG_PYTHON)
  add_subdirectory(python)
endif()

set(WRAP_SWIG_CSHARP OFF CACHE BOOL "Build SWIG CSharp interface")

IF (WRAP_SWIG_CSHARP)
  add_subdirectory(csharp)
endif()



