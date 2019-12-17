set(PROC_PATH ${PROJECT_SOURCE_DIR}/CepGenProcesses)
file(GLOB sources *.cpp Fortran/*.f ${PROC_PATH}/*.cpp)
set(proc_definitions)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(Fortran)

if(ALPHAS_SRC)
  message(STATUS "alphaS evolution algorithm found in ${ALPHAS_SRC}")
  list(APPEND sources ${ALPHAS_SRC})
  list(APPEND proc_definitions -DALPHA_S)
endif()

file(GLOB f77_proc_sources ${PROC_PATH}/*.f ${PROC_PATH}/CepGenWrapper.cpp)
if(f77_proc_sources)
  #message(STATUS "Fortran processes found: ${f77_proc_sources}")
  list(APPEND sources ${f77_proc_sources})
endif()

add_library(CepGenProcesses SHARED ${sources})
target_compile_definitions(CepGenProcesses PUBLIC ${proc_definitions})

install(TARGETS CepGenProcesses DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)
