#--- linking with CTML

set(CTML_DIRS $ENV{CTML_DIR} /usr /usr/local)
find_path(CTML_INCLUDE NAMES CTML/ctml.hpp HINTS ${CTML_DIRS} PATH_SUFFIXES include)

if(NOT CTML_INCLUDE)
  return()
endif()

#----- build the object

cepgen_build(CepGenCTML SOURCES *.cpp
    EXT_HEADERS ${CTML_INCLUDE}
    UTILS utils/*.cc
    INSTALL_COMPONENT doc)
