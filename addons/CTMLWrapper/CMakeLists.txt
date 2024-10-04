#--- linking with CTML

set(CTML_DIRS $ENV{CTML_DIR} /usr /usr/local)
find_path(CTML_INCLUDE NAMES CTML/ctml.hpp HINTS ${CTML_DIRS} PATH_SUFFIXES include)

if(NOT CTML_INCLUDE)
  return()
endif()

#----- build the object
cepgen_build(CepGenCTML SOURCES src/*.cpp
    EXT_HEADERS ${CTML_INCLUDE}
    INSTALL_COMPONENT doc)
