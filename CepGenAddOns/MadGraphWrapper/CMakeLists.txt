#--- linking with ProMC

set(MADGRAPH_DIRS $ENV{MADGRAPH_DIR} /usr /usr/local)
find_program(MADGRAPH_BIN NAMES mg5_aMC HINTS ${MADGRAPH_DIRS} PATH_SUFFIXES bin)

file(GLOB tmpl *.tpp)

#----- build the object

set(mg_defs)
list(APPEND mg_defs "-DMADGRAPH_BIN=\"${MADGRAPH_BIN}\"")
list(APPEND mg_defs "-DMADGRAPH_PROC_TMPL=\"${tmpl}\"")
list(APPEND mg_defs "-DCC_CFLAGS=\"${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -I${PROJECT_SOURCE_DIR}\"")

cepgen_build(CepGenMadGraph
    SOURCES *.cpp
    EXT_LIBS stdc++fs
    DEPENDS ${MADGRAPH_BIN}
    DEFINITIONS ${mg_defs}
    INSTALL_COMPONENT madgraph)
#target_compile_features(CepGenMadGraph PRIVATE cxx_range_for)
