file(GLOB sources *.cpp Fortran/*.f)

include_directories(${PROJECT_SOURCE_DIR})

if(kmr_f_sources)
  list(APPEND sources ${kmr_f_sources})
endif()

add_library(CepGenProcesses SHARED ${sources})

install(TARGETS CepGenProcesses DESTINATION lib)
