#----- searching for Python

set(Python_LIBRARIES)
set(Python_INCLUDE_DIRS)
if(CMAKE_VERSION VERSION_GREATER 3.12)
  find_package(Python COMPONENTS Interpreter Development)
  if(NOT Python_FOUND)
    return()
  endif()
  message(STATUS "Python v${Python_VERSION} found")
else()
  find_package(PythonInterp)
  if(NOT PYTHONINTERP_FOUND)
    return()
  endif()
  find_package(PythonLibs)
  if(NOT PYTHONLIBS_FOUND)
    return()
  endif()
  message(STATUS "Python v${PYTHONLIBS_VERSION_STRING} found")
  list(APPEND ${Python_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
  list(APPEND ${Python_LIBRARIES} ${PYTHON_LIBRARIES})
endif()

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

cepgen_build(CepGenPython SOURCES *.cpp
    EXT_LIBS ${Python_LIBRARIES}
    EXT_HEADERS ${Python_INCLUDE_DIRS}
    INSTALL_COMPONENT python)

