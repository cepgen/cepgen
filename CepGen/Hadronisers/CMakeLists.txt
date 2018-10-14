file(GLOB sources *Hadroniser.cpp PythiaEventInterface.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenExternalHadronisers SHARED ${sources})
target_link_libraries(CepGenExternalHadronisers ${CEPGEN_EXTERNAL_HADR_REQS})

#----- check the external dependencies

install(TARGETS CepGenExternalHadronisers DESTINATION lib)
