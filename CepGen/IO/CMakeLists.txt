file(GLOB io_sources FortranInterface.cpp *Handler.cpp)

include_directories(${PROJECT_SOURCE_DIR})

add_library(CepGenIO SHARED ${io_sources})
target_link_libraries(CepGenIO ${CEPGEN_EXTERNAL_IO_REQS})

install(TARGETS CepGenIO DESTINATION lib)

