file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

add_library(CepGenExporter SHARED ${sources})

install(TARGETS CepGenExporter DESTINATION lib)

