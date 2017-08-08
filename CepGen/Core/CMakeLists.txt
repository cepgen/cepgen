file(GLOB sources *.cpp ../*.cpp)
file(GLOB headers *.h ../*.h)

include_directories(${PROJECT_SOURCE_DIR})

add_library(CepGenCore SHARED ${sources})
target_link_libraries(CepGenCore ${GSL_LIB} ${GSL_CBLAS_LIB} CepGenExternalHadronisers rt)

