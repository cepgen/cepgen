file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

find_package(HepMC)

if(HEPMC_LIBRARIES)
  add_library(Exporter SHARED ${sources})
  target_link_libraries(Exporter ${HEPMC_LIBRARIES})
  add_definitions(-DHEPMC_LINKED)
endif()
