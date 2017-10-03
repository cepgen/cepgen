file(GLOB sources *.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenPhysics SHARED ${sources})

#----- installation rules

install(TARGETS CepGenPhysics DESTINATION lib)

