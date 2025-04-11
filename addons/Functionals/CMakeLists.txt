#--- searching for MuParser
find_library(MUPARSER muparser)
if(MUPARSER)
  cepgen_build(CepGenMuParser
      SOURCES src/FunctionalMuParser.cpp
      LIBRARIES ${MUPARSER}
      COMPONENT muparser)
endif()

#--- searching for exprtk
find_path(EXPRTK exprtk.hpp HINTS $ENV{EXPRTK_DIR} PATH_SUFFIXES include)
if(EXPRTK)
  cepgen_build(CepGenExprtk
      SOURCES src/FunctionalExprTk.cpp
      INCLUDES ${EXPRTK} # header-only
      COMPONENT exprtk)
endif()

#--- searching for TinyExpr
find_path(TINYEXPR_SRC NAMES tinyexpr.c HINTS $ENV{TINYEXPR_DIR})
find_path(TINYEXPR_INCLUDE NAMES tinyexpr.h HINTS $ENV{TINYEXPR_DIR})
if(TINYEXPR_SRC AND TINYEXPR_INCLUDE)
  cepgen_build(CepGenTinyExpr
      SOURCES src/FunctionalTinyExpr.cpp ${TINYEXPR_SRC}/tinyexpr.c
      INCLUDES ${TINYEXPR_INCLUDE}
      COMPONENT tinyexpr)
endif()

#--- searching for libmatheval
find_library(LIBMATHEVAL NAMES matheval HINTS $ENV{LIBMATHEVAL_DIR} PATH_SUFFIXES lib)
find_path(LIBMATHEVAL_INCLUDE NAMES matheval.h HINTS $ENV{LIBMATHEVAL_DIR} PATH_SUFFIXES include)
if(LIBMATHEVAL)
  cepgen_build(CepGenLibmatheval
      SOURCES src/FunctionalLibmatheval.cpp
      LIBRARIES ${LIBMATHEVAL}
      INCLUDES ${LIBMATHEVAL_INCLUDE}
      COMPONENT libmatheval)
endif()

#--- searching for ATMSP
find_path(ATMSP_INCLUDE NAMES atmsp.h HINTS $ENV{ATMSP_DIR} PATH_SUFFIXES include)
if(ATMSP_INCLUDE)
  cepgen_build(CepGenATMSP
      SOURCES src/FunctionalATMSP.cpp
      INCLUDES ${ATMSP_INCLUDE}
      COMPONENT atmsp)
endif()
