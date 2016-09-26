file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

find_package(HepMC REQUIRED)

add_library(Exporter OBJECT ${sources})
