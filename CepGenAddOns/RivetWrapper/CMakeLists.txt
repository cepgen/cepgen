#--- searching for YODA
find_library(YODA YODA HINTS $ENV{YODA_DIR} PATH_SUFFIXES lib)
find_path(YODA_INCLUDE YODA HINTS $ENV{YODA_DIR} PATH_SUFFIXES include)

if(NOT YODA)
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources YODAHistsHandler.cpp)

#----- build the object

cepgen_build(CepGenRivet SOURCES YODAHistsHandler.cpp
    EXT_LIBS ${YODA}
    EXT_HEADERS ${YODA_INCLUDE}
    INSTALL_COMPONENT rivet)
