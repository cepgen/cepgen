#----- searching for Python

set(Python_LIBRARIES)
set(Python_INCLUDE_DIRS)
find_package(Python COMPONENTS Interpreter Development)
if(NOT Python_FOUND)
  message(WARNING "Python interpreter+development packages were not found")
  return()
endif()

message(NOTICE "Python v${Python_VERSION} found")
if(${Python_VERSION_MAJOR} VERSION_GREATER_EQUAL 3)
  set(tests test/*.cc)
endif()

find_package(pybind11 CONFIG QUIET)
if(pybind11_FOUND)  # build the pybind11 Python binding module
  message(NOTICE "Python bindings will be built")
  pybind11_add_module(pycepgen NO_EXTRAS src/PythonWrapper.cpp src/PythonUtils.cpp)
  target_link_libraries(pycepgen PRIVATE CepGen)
  set_target_properties(pycepgen PROPERTIES CXX_VISIBILITY_PRESET default VISIBILITY_INLINES_HIDDEN OFF)
  target_include_directories(pycepgen PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  install(TARGETS pycepgen DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT python)
endif()

#----- build the object
enable_testing()
cepgen_build(CepGenPython
    SOURCES src/CardHandler.cpp
            src/ConfigWriter.cpp
            src/DocumentationGenerator.cpp
            src/Environment.cpp
            src/Error.cpp
            src/Functional.cpp
            src/Integrator.cpp
            src/ObjectPtr.cpp
            src/Utils.cpp
    TESTS ${tests}
    LIBRARIES ${Python_LIBRARIES}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${Python_INCLUDE_DIRS}
    COMPONENT python)
install(DIRECTORY CepGenPython
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT python)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/python
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/CepGen
    COMPONENT python
    FILES_MATCHING PATTERN "*.py"
    PATTERN "*.pyc" EXCLUDE)

# Python description of modules
if(CMAKE_BUILD_UTILS)
  set(MODULES_DESCRIPTOR "${PROJECT_BINARY_DIR}/bin/cepgenDescribeModules")
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

# Python tests
if(CMAKE_BUILD_TESTS)
  add_test(NAME pycepgen_module
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/pycepgen_module.py)
  set_property(TEST pycepgen_module PROPERTY ENVIRONMENT
    "CEPGEN_PATH=${PROJECT_SOURCE_DIR};PYTHONPATH=$<TARGET_FILE_DIR:pycepgen>")
endif()
