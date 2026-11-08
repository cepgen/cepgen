#--- searching for LHAPDF
find_library(LHAPDF LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES lib)
find_path(LHAPDF_INCLUDE LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES include)

if(NOT LHAPDF)
  return()
endif()

#----- build the object
cepgen_build(CepGenLHAPDF
    SOURCES src/*.cpp
    LIBRARIES ${LHAPDF}
    INCLUDES ${LHAPDF_INCLUDE}
    COMPONENT lhapdf)
cpack_add_component(lhapdf
    DISPLAY_NAME "CepGen LHAPDF wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the LHAPDF library"
    DEPENDS lib)
