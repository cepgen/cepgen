#--- linking with ProMC

set(MADGRAPH_DIRS $ENV{MADGRAPH_DIR} /usr /usr/local)
find_program(MADGRAPH_BIN NAMES mg5_aMC HINTS ${MADGRAPH_DIRS} PATH_SUFFIXES bin)

if(MADGRAPH_BIN)
  message(STATUS "MadGraph_aMC executable found in ${MADGRAPH_BIN}")
else()
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources *.cpp)
file(GLOB tmpl *.tpp)

#----- build the objects

set(mg_defs)
list(APPEND mg_defs "-DMADGRAPH_BIN=\"${MADGRAPH_BIN}\"")
list(APPEND mg_defs "-DMADGRAPH_PROC_TMPL=\"${tmpl}\"")
list(APPEND mg_defs "-DCC_CFLAGS=\"${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -I${PROJECT_SOURCE_DIR}\"")
message(STATUS ${mg_defs})

add_library(CepGenMadGraph SHARED ${sources})
#target_compile_features(CepGenMadGraph PRIVATE cxx_range_for)
target_link_libraries(CepGenMadGraph stdc++fs)
target_compile_definitions(CepGenMadGraph PRIVATE ${mg_defs})

#----- installation rules

install(TARGETS CepGenMadGraph DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT madgraph)
