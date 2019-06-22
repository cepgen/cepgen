file(GLOB core_sources Core/*.cpp *.cpp)
file(GLOB phys_sources Physics/*.cpp)
file(GLOB sf_sources StructureFunctions/*.cpp)
file(GLOB io_sources IO/*.cpp)
file(GLOB hadr_sources Hadronisers/*Hadroniser.cpp Hadronisers/PythiaEventInterface.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- check the external dependencies for SFs

set(GRV_PATH ${PROJECT_SOURCE_DIR}/External)
file(GLOB grv_sources ${GRV_PATH}/grv_*.f)
if(grv_sources)
  message(STATUS "GRV PDFset found in ${grv_sources}!")
  add_definitions(-DGRVPDF)
  list(APPEND sf_sources ${grv_sources})
else()
  message(STATUS "GRV PDFset not found. Will proceed without it")
endif()

#----- check the external dependencies for physics utilities

if(alphas_sources)
  list(APPEND phys_sources ${alphas_sources})
endif()

#--- searching for Pythia 6
if(EXISTS $ENV{PYTHIA6_SRC})
  message(STATUS "Pythia 6 source found in $ENV{PYTHIA6_SRC}")
  file(GLOB pythia6_src $ENV{PYTHIA6_SRC})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-tabs -Wno-maybe-uninitialized -Wno-integer-division -Wno-unused-variable -Wno-unused-dummy-argument")
  add_library(Pythia6 SHARED ${pythia6_src})
  list(APPEND CEPGEN_EXTERNAL_HADR_REQS Pythia6)
endif()

#----- build the objects

add_library(CepGenCore SHARED ${core_sources} ${phys_sources} ${sf_sources})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_CORE_REQS})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_STRF_REQS})
add_library(CepGenAddOns SHARED ${io_sources} ${hadr_sources})
target_link_libraries(CepGenAddOns ${CEPGEN_EXTERNAL_IO_REQS})
target_link_libraries(CepGenAddOns ${CEPGEN_EXTERNAL_HADR_REQS})

#----- installation rules

install(TARGETS CepGenCore DESTINATION lib)
install(TARGETS CepGenAddOns DESTINATION lib)

