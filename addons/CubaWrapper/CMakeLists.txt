#--- searching for Cuba
find_path(CUBA_INCLUDE tools HINTS $ENV{CUBA_DIR})
find_library(CUBA cuba HINTS $ENV{CUBA_DIR})
# Cuba must be compiled with fPIC, i.e.
# ./configure CFLAGS=-fPIC CXXFLAGS=-fPIC

if(NOT CUBA)
  return()
endif()

#----- build the object
cepgen_build(CepGenCuba SOURCES src/*.cpp
    EXT_LIBS ${CUBA}
    EXT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR} ${CUBA_INCLUDE}
    PROPERTY POSITION_INDEPENDENT_CODE ON
    INSTALL_DIRECTORY CepGenCuba
    INSTALL_COMPONENT cuba)

#set_target_properties(CepGenCore PROPERTIES COMPILE_FLAGS "-ffast-math -lm")
