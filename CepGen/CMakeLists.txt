file(GLOB core_sources Core/*.cpp *.cpp)
file(GLOB cards_sources Cards/*.cpp)

include_directories(${PROJECT_SOURCE_DIR})

add_library(CepGenCore SHARED ${core_sources} ${cards_sources})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_CORE_REQS})

install(TARGETS CepGenCore DESTINATION lib)

