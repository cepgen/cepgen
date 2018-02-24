file(GLOB sources *.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenExternalHadronisers SHARED ${sources})
if(PYTHIA8)
  target_link_libraries(CepGenExternalHadronisers ${PYTHIA8} dl)
endif()

#----- check the external dependencies

install(TARGETS CepGenExternalHadronisers DESTINATION lib)
