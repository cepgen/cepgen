file(GLOB core_sources *.cpp ../*.cpp)
file(GLOB cards_sources ${CEPGEN_SOURCE_DIR}/Cards/*.cpp)

include_directories(${PROJECT_SOURCE_DIR})

add_library(CepGenCore SHARED ${core_sources} ${cards_sources})
target_link_libraries(CepGenCore rt)

install(TARGETS CepGenCore DESTINATION lib)

