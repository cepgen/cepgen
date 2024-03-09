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
