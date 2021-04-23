include_directories(${PROJECT_SOURCE_DIR})

file(GLOB proc_sources *.cpp *.f)

#----- build the objects

add_library(CepGenProcesses SHARED ${proc_sources})

#----- installation rules

install(TARGETS CepGenProcesses DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lib)
