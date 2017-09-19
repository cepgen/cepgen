file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenExternalHadronisers SHARED ${sources})

#----- check the external dependencies

install(TARGETS CepGenExternalHadronisers DESTINATION lib)
