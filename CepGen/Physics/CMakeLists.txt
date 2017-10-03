file(GLOB phys_sources *.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenPhysics SHARED ${phys_sources})

#----- installation rules

install(TARGETS CepGenPhysics DESTINATION lib)

