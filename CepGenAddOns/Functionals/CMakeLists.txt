set(CEPGEN_FUNCTIONALS_SRC)
set(CEPGEN_FUNCTIONALS_LIBS)
set(CEPGEN_FUNCTIONALS_INCLUDES)

find_library(MUPARSER muparser)
if(MUPARSER)
  message(STATUS "Formulas parsing will use muParser found in ${MUPARSER}")
  list(APPEND CEPGEN_FUNCTIONALS_LIBS ${MUPARSER})
  file(GLOB func_muparser FunctionalMuParser.cpp)
  list(APPEND CEPGEN_FUNCTIONALS_SRC ${func_muparser})
endif()

find_path(EXPRTK exprtk.hpp HINTS $ENV{EXPRTK_DIR} PATH_SUFFIXES include)
if(EXPRTK)
  message(STATUS "Formulas parsing will use exprtk found in ${EXPRTK}")
  list(APPEND CEPGEN_FUNCTIONALS_INCLUDES ${EXPRTK}) # header-only
  file(GLOB func_exprtk FunctionalExprTk.cpp)
  list(APPEND CEPGEN_FUNCTIONALS_SRC ${func_exprtk})
endif()

if(NOT CEPGEN_FUNCTIONALS_SRC)
  return()
endif()

cepgen_build(CepGenFunctionals
        SOURCES ${CEPGEN_FUNCTIONALS_SRC}
        EXT_LIBS ${CEPGEN_FUNCTIONALS_LIBS}
        EXT_HEADERS ${CEPGEN_FUNCTIONALS_INCLUDES}
        INSTALL_COMPONENT lib)
