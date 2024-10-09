set(CLHEP_DIRS $ENV{CLHEP_DIR} /usr /usr/local)
find_package(CLHEP HINTS ${CLHEP_DIRS})
if(NOT CLHEP_FOUND)
  return()
endif()

#----- build the object
cepgen_build(CepGenCLHEP SOURCES src/*.cpp
    LIBRARIES CLHEP::CLHEP
    COMPONENT clhep)
cpack_add_component(clhep
    DISPLAY_NAME "CepGen CLHEP wrappers library"
    DESCRIPTION "A CepGen wrapper to the CLHEP library"
    DEPENDS lib)
