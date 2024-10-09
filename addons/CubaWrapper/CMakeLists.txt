#--- searching for Cuba
find_path(CUBA_INCLUDE tools HINTS $ENV{CUBA_DIR})
find_library(CUBA cuba HINTS $ENV{CUBA_DIR})
# Cuba must be compiled with fPIC, i.e.
# ./configure CFLAGS=-fPIC CXXFLAGS=-fPIC

if(NOT CUBA)
  return()
endif()

#----- build the object
cepgen_build(CepGenCuba
    SOURCES src/*.cpp
    LIBRARIES ${CUBA}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${CUBA_INCLUDE}
    PROPERTY POSITION_INDEPENDENT_CODE ON
    COMPONENT cuba)
install(DIRECTORY CepGenCuba
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT cuba)

#set_target_properties(CepGenCore PROPERTIES COMPILE_FLAGS "-ffast-math -lm")
