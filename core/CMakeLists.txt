file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

find_package(GSL REQUIRED)

add_library(CepGenCore SHARED ${sources})
target_link_libraries(CepGenCore ${GSL_LIBRARIES} rt)

