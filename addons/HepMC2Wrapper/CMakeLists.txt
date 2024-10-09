#--- searching for HepMC2
set(HEPMC_DIRS $ENV{HEPMC_DIR} /usr /usr/local)
find_library(HEPMC_LIB NAME HepMC HINTS ${HEPMC_DIRS} PATH_SUFFIXES lib64 lib)
find_path(HEPMC_INCLUDE NAME HepMC HINTS ${HEPMC_DIRS} PATH_SUFFIXES include)

if(NOT HEPMC_LIB)
  return()
endif()

#----- build the object
cepgen_build(CepGenHepMC2
    SOURCES src/*.cpp
    LIBRARIES ${HEPMC_LIB}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${HEPMC_INCLUDE}
    TESTS test/*.cc
    COMPONENT hepmc2)
install(DIRECTORY CepGenHepMC2
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT hepmc2)
cpack_add_component(hepmc2
    DISPLAY_NAME "CepGen HepMC wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the HepMC library"
    DEPENDS lib)
