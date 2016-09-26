file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

find_package(HepMC REQUIRED)

if(${HEPMC_FOUND})
  add_library(Exporter SHARED ${sources})
  target_link_libraries(Exporter ${HEPMC_LIBRARIES})
endif()
