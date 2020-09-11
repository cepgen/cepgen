#--- searching for YODA
find_library(YODA YODA HINTS $ENV{YODA_DIR} PATH_SUFFIXES lib)
find_path(YODA_INCLUDE YODA HINTS $ENV{YODA_DIR} PATH_SUFFIXES include)

if(NOT YODA)
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources YODAHistsHandler.cpp)

#----- build the objects

add_library(CepGenRivet SHARED ${sources})
target_link_libraries(CepGenRivet PRIVATE ${YODA})
target_include_directories(CepGenRivet PRIVATE ${YODA_INCLUDE})
list(APPEND CEPGEN_ADDONS_EXT CepGenRivet)
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenRivet DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT rivet)
