include_directories(${PROJECT_SOURCE_DIR}/CepGen/Process/Fortran)

#--- build the processes matrix elements
file(GLOB proc_sources *.cpp *.f)
cepgen_build(CepGenProcesses
        SOURCES ${proc_sources}
        INSTALL_COMPONENT lib)

#--- also build examples if requested
if(CMAKE_BUILD_PROCESSES_EXAMPLES)
  file(GLOB ex_proc_sources Examples/*.cpp Examples/*.f)
  cepgen_build(CepGenProcessesExamples
               SOURCES ${ex_proc_sources})
endif()
