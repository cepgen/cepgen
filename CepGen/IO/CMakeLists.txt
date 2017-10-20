file(GLOB io_sources FortranInterface.cpp)

include_directories(${PROJECT_SOURCE_DIR})

add_library(CepGenIO SHARED ${io_sources})

install(TARGETS CepGenIO DESTINATION lib)

