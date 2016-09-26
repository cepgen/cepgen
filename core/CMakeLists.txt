file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

find_package(GSL REQUIRED)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -gfortran")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -lgsl -lgslcblas")

add_library(Core SHARED ${sources})
