file(GLOB ev_sources *.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenEvent SHARED ${ev_sources})

#----- installation recipes

install(TARGETS CepGenEvent DESTINATION ${LIB_INSTALL_DIR} COMPONENT lib)

