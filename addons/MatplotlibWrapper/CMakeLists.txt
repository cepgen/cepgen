#--- searching for MATPLOTLIB
find_package(matplotlib_cpp QUIET)

if(NOT ${matplotlib_cpp_FOUND})
  return()
endif()

set(Python_LIBRARIES)
set(Python_INCLUDE_DIRS)
if(CMAKE_VERSION VERSION_GREATER 3.12)
  find_package(Python COMPONENTS Interpreter Development)
  if(NOT Python_FOUND)
    return()
  endif()
else()
  find_package(PythonInterp)
  if(NOT PYTHONINTERP_FOUND)
    return()
  endif()
  find_package(PythonLibs)
  if(NOT PYTHONLIBS_FOUND)
    return()
  endif()
  list(APPEND Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
  list(APPEND Python_LIBRARIES ${PYTHON_LIBRARIES})
endif()

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object
cepgen_build(CepGenMatplotlib SOURCES src/*.cpp
    EXT_LIBS ${matplotlib_cpp_LIBRARIES} ${Python_LIBRARIES}
    EXT_HEADERS ${matplotlib_cpp_INCLUDE_DIRS} ${Python_INCLUDE_DIRS}
    INSTALL_COMPONENT matplotlib)
