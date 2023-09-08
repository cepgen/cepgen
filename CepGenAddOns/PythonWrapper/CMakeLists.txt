#----- searching for Python

set(Python_LIBRARIES)
set(Python_INCLUDE_DIRS)
if(CMAKE_VERSION VERSION_GREATER 3.12)
  find_package(Python COMPONENTS Interpreter Development)
  if(NOT Python_FOUND)
    return()
  endif()
  message(STATUS "Python v${Python_VERSION} found")
  if(${Python_VERSION_MAJOR} VERSION_GREATER_EQUAL 3)
    file(GLOB tests test/*.cc)
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
  message(STATUS "Python v${PYTHONLIBS_VERSION_STRING} found")
  list(APPEND Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
  list(APPEND Python_LIBRARIES ${PYTHON_LIBRARIES})
endif()

file(GLOB sources *.cpp)

execute_process(COMMAND pip show vegas RESULT_VARIABLE VEGAS_RES OUTPUT_QUIET)
if(${VEGAS_RES} EQUAL 0)
  message(STATUS "... Vegas python package present")
else()
  file(GLOB vegas IntegratorVegasPlus.cpp)
  list(REMOVE_ITEM sources ${vegas})
endif()

#----- build the object

enable_testing()
cepgen_build(CepGenPython SOURCES ${sources}
    TESTS ${tests}
    EXT_LIBS ${Python_LIBRARIES}
    EXT_HEADERS ${Python_INCLUDE_DIRS}
    INSTALL_COMPONENT python)
