set(GNUPLOT_DIRS $ENV{GNUPLOT_DIR} /usr /usr/local)
find_program(GNUPLOT_BIN NAMES gnuplot HINTS ${GNUPLOT_DIRS} PATH_SUFFIXES bin)
if(NOT GNUPLOT_BIN)
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object
cepgen_build(CepGenGnuplot SOURCES src/*.cpp
    DEPENDS ${GNUPLOT_BIN}
    DEFINITIONS "-DGNUPLOT_BIN=${GNUPLOT_BIN}"
    COMPONENT gnuplot)
cpack_add_component(gnuplot
    DISPLAY_NAME "CepGen Gnuplot wrappers library"
    DESCRIPTION "A CepGen wrapper to the Gnuplot plotter"
    DEPENDS lib)
