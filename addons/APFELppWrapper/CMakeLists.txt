#--- searching for APFEL++
set(apfelxx_DIRS $ENV{apfelxx_DIR} /usr /usr/local)
find_library(apfelxx apfelxx HINTS ${apfelxx_DIRS} PATH_SUFFIXES lib)
find_path(apfelxx_INCLUDE apfel/apfelxx.h HINTS ${apfelxx_DIRS} PATH_SUFFIXES include)

if(NOT apfelxx)
  return()
endif()

#----- build the object
cepgen_build(CepGenAPFELpp
    SOURCES src/*.cpp
    LIBRARIES ${apfelxx}
    INCLUDES ${apfelxx_INCLUDE}
    COMPONENT apfelxx)
cpack_add_component(apfelxx
    DISPLAY_NAME "CepGen APFEL++ wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the APFEL++ library"
    DEPENDS lib)
