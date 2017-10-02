file(GLOB core_sources Core/*.cpp *.cpp)
file(GLOB cards_sources Cards/*.cpp)
file(GLOB ev_sources Event/*.cpp)

include_directories(${PROJECT_SOURCE_DIR})

add_library(CepGenCore SHARED ${core_sources} ${cards_sources} ${ev_sources})
target_link_libraries(CepGenCore rt)

install(TARGETS CepGenCore DESTINATION lib)

