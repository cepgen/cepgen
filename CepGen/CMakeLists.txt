file(GLOB core_sources Core/*.cpp *.cpp)
file(GLOB phys_sources Physics/*.cpp)
file(GLOB sf_sources StructureFunctions/*.cpp)
file(GLOB io_sources IO/*.cpp)
file(GLOB hadr_sources Hadronisers/GenericHadroniser.cpp)

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

set(hadr_libraries "")

#--- linking with Pythia 8

if(PYTHIA8)
  message(STATUS "Pythia 8 found in ${PYTHIA8}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
  include_directories(${PYTHIA8_INCLUDE})
  add_definitions(-DPYTHIA8)
  list(APPEND hadr_sources "Hadronisers/Pythia8Hadroniser.cpp")
  list(APPEND hadr_sources "Hadronisers/PythiaEventInterface.cpp")
  list(APPEND hadr_libraries ${PYTHIA8} dl)
endif()

#----- build the objects

add_library(CepGenCore SHARED ${core_sources} ${phys_sources} ${sf_sources})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_CORE_REQS})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_STRF_REQS})

add_library(CepGenAddOns SHARED ${io_sources} ${hadr_sources})
target_link_libraries(CepGenAddOns ${hadr_libraries})
target_link_libraries(CepGenAddOns ${CEPGEN_EXTERNAL_IO_REQS})
target_link_libraries(CepGenAddOns CepGenEvent)

#----- installation rules

install(TARGETS CepGenCore DESTINATION lib)
install(TARGETS CepGenAddOns DESTINATION lib)

