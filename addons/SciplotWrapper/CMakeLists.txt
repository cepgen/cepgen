find_package(sciplot)
if(NOT sciplot_FOUND)
  return()
endif()

cepgen_build(CepGenSciplot SOURCES src/*.cpp
    LIBRARIES sciplot::sciplot
    COMPONENT sciplot)
cpack_add_component(sciplot
    DISPLAY_NAME "CepGen Sciplot wrappers library"
    DESCRIPTION "A CepGen wrapper to the Sciplot plotter"
    DEPENDS lib)
