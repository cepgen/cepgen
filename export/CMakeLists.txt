file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

find_package(HepMC)

if(HEPMC_LIBRARIES)
  add_library(CepGenExporter SHARED ${sources})
  target_link_libraries(CepGenExporter ${HEPMC_LIBRARIES})
  add_definitions(-DHEPMC_LINKED)
endif()
