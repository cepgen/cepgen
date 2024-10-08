#--- searching for Pythia 6
set(PYTHIA6_DIRS $ENV{PYTHIA6_DIR} /usr /usr/local /opt/pythia6)
find_library(PYTHIA6 pythia6 HINTS ${PYTHIA6_DIRS} PATH_SUFFIXES lib)

set(PYTHIA6LIB)
if(PYTHIA6)
  list(APPEND PYTHIA6LIB ${PYTHIA6})
  find_library(PYTHIA6DUMMY pythia6_dummy HINTS ${PYTHIA6_DIRS} PATH_SUFFIXES lib)
  if(PYTHIA6DUMMY)
    list(APPEND PYTHIA6LIB ${PYTHIA6DUMMY})
  endif()
  find_library(PYTHIA6PDFDUMMY pythia6_pdfdummy HINTS ${PYTHIA6_DIRS} PATH_SUFFIXES lib)
  if(PYTHIA6PDFDUMMY)
    list(APPEND PYTHIA6LIB ${PYTHIA6PDFDUMMY})
  endif()
else()
  file(GLOB PYTHIA6_SRC ${PROJECT_SOURCE_DIR}/External/pythia*.f)
  if(EXISTS $ENV{PYTHIA6_SRC})
    file(GLOB pythia6_src $ENV{PYTHIA6_SRC})
  elseif(EXISTS ${PYTHIA6_SRC})
    file(GLOB pythia6_src ${PYTHIA6_SRC})
  else()
    return()
  endif()
  message(STATUS "Pythia 6 source found in ${pythia6_src}")
  add_library(pythia6 SHARED ${pythia6_src})
  install(TARGETS pythia6 DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT pythia6)
  target_compile_options(pythia6 PRIVATE "-Wno-tabs -Wno-maybe-uninitialized -Wno-integer-division -Wno-unused-variable -Wno-unused-dummy-argument")
  list(APPEND PYTHIA6LIB pythia6)
endif()

enable_testing()

#----- build the objects
cepgen_build(CepGenPythia6
    SOURCES src/*.cpp
    TESTS test/*.cc
    LIBRARIES ${PYTHIA6LIB}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    OPTIONS "-Wno-misleading-indentation"
    COMPONENT pythia6)
install(DIRECTORY CepGenPythia6
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT pythia6)
add_test(NAME pythia6_kt_hadr COMMAND ${CMAKE_BINARY_DIR}/test/generator/ktfactor_proc_hadronisation -H pythia6)
