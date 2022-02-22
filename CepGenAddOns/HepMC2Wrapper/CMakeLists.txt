#--- searching for HepMC2
set(HEPMC_DIRS $ENV{HEPMC_DIR} /usr /usr/local)
find_library(HEPMC_LIB NAME HepMC HINTS ${HEPMC_DIRS} PATH_SUFFIXES lib64 lib)
find_path(HEPMC_INCLUDE NAME HepMC HINTS ${HEPMC_DIRS} PATH_SUFFIXES include)

if(NOT HEPMC_LIB)
  return()
endif()

#----- build the object

cepgen_build(CepGenHepMC2 SOURCES *.cpp
    EXT_LIBS ${HEPMC_LIB}
    EXT_HEADERS ${HEPMC_INCLUDE}
    INSTALL_COMPONENT hepmc2)
