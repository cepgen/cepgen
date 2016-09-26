file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

find_package(GSL REQUIRED)

add_library(Core SHARED ${sources})
target_link_libraries(Core ${GSL_LIBRARIES} rt)
