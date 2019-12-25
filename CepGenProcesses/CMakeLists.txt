include_directories(${PROJECT_SOURCE_DIR})

file(GLOB proc_sources *.cpp *.f)

#----- build the objects

add_library(CepGenProcesses SHARED ${proc_sources})
target_include_directories(CepGenProcesses PUBLIC ${CEPGEN_LIB_DIR}/Processes/Fortran)

#----- installation rules

install(TARGETS CepGenProcesses DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)
