set(MADGRAPH_DIRS $ENV{MADGRAPH_DIR} /usr /usr/local)
find_program(MADGRAPH_BIN NAMES mg5_aMC HINTS ${MADGRAPH_DIRS} PATH_SUFFIXES bin)
if(NOT MADGRAPH_BIN)
  return()
endif()

#----- build the object

set(mg_defs)
file(GLOB tmpl template/*.cpp)
list(APPEND mg_defs "-DMADGRAPH_BIN=\"${MADGRAPH_BIN}\"")
list(APPEND mg_defs "-DMADGRAPH_PROC_TMPL=\"${tmpl}\"")
list(APPEND mg_defs "-DCC_CFLAGS=\"${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -I${PROJECT_SOURCE_DIR}\"")

cepgen_build(CepGenMadGraph
    SOURCES *.cpp
    EXT_LIBS stdc++fs CepGenPython
    EXT_HEADERS ${Python_INCLUDE_DIRS}
    DEPENDS ${MADGRAPH_BIN}
    DEFINITIONS ${mg_defs}
    TESTS test/*.cc
    INSTALL_COMPONENT madgraph)

add_test(NAME cepgen_mg5 COMMAND cepgen Cards/mg5amcnlo/aamumu_cfg.py -n 1000)
set_property(TEST cepgen_mg5 PROPERTY ENVIRONMENT "CEPGEN_PATH=${CEPGEN_PATH}")
