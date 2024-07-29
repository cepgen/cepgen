include_directories(${PROJECT_SOURCE_DIR}/CepGen/Process/Fortran)

#--- build the processes matrix elements
file(GLOB proc_sources *.cpp *.f)
cepgen_build(CepGenProcesses
        SOURCES ${proc_sources}
        INSTALL_COMPONENT lib)
