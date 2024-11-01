set(asciichart_DIRS $ENV{asciichart_DIR} /usr /usr/local)
find_path(asciichart_INCLUDE_DIR NAMES ascii/ascii.h HINTS ${asciichart_DIRS} PATH_SUFFIXES include)
if(NOT asciichart_INCLUDE_DIR)
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object
cepgen_build(CepGenasciiChart SOURCES src/*.cpp
    INCLUDES ${asciichart_INCLUDE_DIR}
    COMPONENT asciichart)
cpack_add_component(asciichart
    DISPLAY_NAME "CepGen asciichart wrappers library"
    DESCRIPTION "A CepGen wrapper to the asciichart plotter"
    DEPENDS lib)
