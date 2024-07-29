#----- searching for Python

set(Python_LIBRARIES)
set(Python_INCLUDE_DIRS)
if(CMAKE_VERSION VERSION_GREATER 3.12)
  find_package(Python COMPONENTS Interpreter Development)
  if(NOT Python_FOUND)
    message(WARNING "Python interpreter+development packages were not found")
    return()
  endif()
  message(NOTICE "Python v${Python_VERSION} found")
  if(${Python_VERSION_MAJOR} VERSION_GREATER_EQUAL 3)
    set(tests test/*.cc)
  endif()
else()
  find_package(PythonInterp)
  if(NOT PYTHONINTERP_FOUND)
    message(WARNING "Python interpreter was not found")
    return()
  endif()
  find_package(PythonLibs)
  if(NOT PYTHONLIBS_FOUND)
    message(WARNING "Python libraries were not found")
    return()
  endif()
  message(STATUS "Python v${PYTHONLIBS_VERSION_STRING} found")
  list(APPEND Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
  list(APPEND Python_LIBRARIES ${PYTHON_LIBRARIES})
endif()

#----- build the object

enable_testing()
cepgen_build(CepGenPython SOURCES *.cpp
    TESTS ${tests}
    EXT_LIBS ${Python_LIBRARIES}
    EXT_HEADERS ${Python_INCLUDE_DIRS}
    INSTALL_COMPONENT python)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/python
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/CepGen
  COMPONENT python
  FILES_MATCHING PATTERN "*.py"
  PATTERN "*.pyc" EXCLUDE)
if(CMAKE_BUILD_UTILS)
  set(MODULES_DESCRIPTOR ${PROJECT_BINARY_DIR}/bin/cepgenDescribeModules)
  set(MODULES_OUTPUT ${PROJECT_BINARY_DIR}/python_modules)
  add_custom_target(CepGenPythonModules ALL
    DEPENDS ${MODULES_DESCRIPTOR}
    COMMENT "Generation of Python modules descriptions")
  add_custom_command(TARGET CepGenPythonModules POST_BUILD
    COMMAND ${CMAKE_COMMAND}
      -DMODULE_DOC_EXE=${MODULES_DESCRIPTOR}
      -DMODULE_DOC_OUTPUT=${MODULES_OUTPUT}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GenerateModulesDocs.cmake)
  install(DIRECTORY ${MODULES_OUTPUT}
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/CepGen
    COMPONENT python
    FILES_MATCHING PATTERN "*.py"
    PATTERN "__pycache__" EXCLUDE)
endif()

cepgen_test_declare_unstable(python_output python_particle_mod)
