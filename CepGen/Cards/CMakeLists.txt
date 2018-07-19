file(GLOB cards_sources *.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenCards SHARED ${cards_sources})
target_link_libraries(CepGenCards ${CEPGEN_EXTERNAL_CARDS_REQS})

#----- installation rules

install(TARGETS CepGenCards DESTINATION lib)

