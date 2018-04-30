file(GLOB sources *.cpp Fortran/*.f)

include_directories(${PROJECT_SOURCE_DIR})

add_library(CepGenProcesses SHARED ${sources})

install(TARGETS CepGenProcesses DESTINATION lib)
