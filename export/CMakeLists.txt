file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

find_package(HepMC REQUIRED)

if(NOT DEFINED ${HEPMC_LIB})
  add_library(Exporter SHARED ${sources})
  target_link_libraries(Exporter ${HEPMC_LIB})
  add_definitions(-DHEPMC_LINKED)
endif()
