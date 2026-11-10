#--- searching for NeoPDF
find_library(NeoPDF neopdf_capi HINTS $ENV{NeoPDF_PATH} PATH_SUFFIXES lib64 lib)
find_path(NeoPDF_INCLUDE NeoPDF.hpp HINTS $ENV{NeoPDF_PATH} PATH_SUFFIXES include/neopdf_capi include)

if(NOT NeoPDF OR NOT NeoPDF_INCLUDE)
  return()
endif()

#----- build the object
cepgen_build(CepGenNeoPDF
    SOURCES src/*.cpp
    LIBRARIES ${NeoPDF}
    INCLUDES ${NeoPDF_INCLUDE}
    COMPONENT neopdf)
cpack_add_component(neopdf
    DISPLAY_NAME "CepGen NeoPDF wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the NeoPDF library"
    DEPENDS lib)
