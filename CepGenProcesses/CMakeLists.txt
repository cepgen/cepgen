include_directories(${PROJECT_SOURCE_DIR})

file(GLOB proc_sources *.cpp *.f)
cepgen_build(CepGenProcesses
        SOURCES ${proc_sources}
        INSTALL_COMPONENT lib)
