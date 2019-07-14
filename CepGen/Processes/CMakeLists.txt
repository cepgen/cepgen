file(GLOB sources *.cpp Fortran/*.f)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(Fortran)

set(F77_PROC_PATH ${PROJECT_SOURCE_DIR}/External/Processes)
file(GLOB oth_proc_sources ${F77_PROC_PATH}/*.f ${F77_PROC_PATH}/CepGenWrapper.cpp)
if(oth_proc_sources)
  message(STATUS "Fortran processes found: ${oth_proc_sources}")
  list(APPEND sources ${oth_proc_sources})
endif()

add_library(CepGenProcesses SHARED ${sources})

install(TARGETS CepGenProcesses DESTINATION lib)
