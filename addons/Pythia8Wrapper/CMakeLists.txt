#--- searching for Pythia 8
set(PYTHIA8_DIRS $ENV{PYTHIA8_DIR} /usr /usr/local /opt/pythia8)
find_library(PYTHIA8 pythia8 HINTS ${PYTHIA8_DIRS} PATH_SUFFIXES lib lib64)
find_path(PYTHIA8_INCLUDE Pythia8 HINTS ${PYTHIA8_DIRS} PATH_SUFFIXES include/Pythia8 include/pythia8 include)

if(NOT PYTHIA8)
  return()
endif()

#--- searching for gzip
set(PYTHIA_DEFS)
find_program(GZIP_BIN gzip DOC "Path to gzip program")
if(GZIP_BIN)
  list(APPEND PYTHIA_DEFS "-DGZIP_BIN=\"${GZIP_BIN}\"")
endif()

enable_testing()

#----- build the object
cepgen_build(CepGenPythia8
    SOURCES src/*.cpp
    LIBRARIES ${PYTHIA8} dl
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHIA8_INCLUDE}
    DEFINITIONS ${PYTHIA_DEFS}
    TESTS test/*.cc
    OPTIONS "-Wno-misleading-indentation"
    COMPONENT pythia8)
cpack_add_component(pythia8
    DISPLAY_NAME "CepGen Pythia 8 wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to Pythia 8"
    DEPENDS lib)
add_test(NAME pythia8_kt_hadr COMMAND ${CMAKE_BINARY_DIR}/test/generator/ktfactor_proc_hadronisation -H pythia8)
