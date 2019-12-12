set(PROC_PATH ${PROJECT_SOURCE_DIR}/CepGenProcesses)
file(GLOB sources *.cpp Fortran/*.f ${PROC_PATH}/*.cpp)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(Fortran)

file(GLOB f77_proc_sources ${PROC_PATH}/*.f ${PROC_PATH}/CepGenWrapper.cpp)
if(f77_proc_sources)
  #message(STATUS "Fortran processes found: ${f77_proc_sources}")
  list(APPEND sources ${f77_proc_sources})
endif()

add_library(CepGenProcesses SHARED ${sources})

install(TARGETS CepGenProcesses DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)
