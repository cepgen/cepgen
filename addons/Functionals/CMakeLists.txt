#--- searching for MuParser
find_library(MUPARSER muparser)
if(MUPARSER)
  cepgen_build(CepGenMuParser
      SOURCES src/FunctionalMuParser.cpp
      EXT_LIBS ${MUPARSER}
      INSTALL_COMPONENT muparser)
endif()

#--- searching for exprtk
find_path(EXPRTK exprtk.hpp HINTS $ENV{EXPRTK_DIR} PATH_SUFFIXES include)
if(EXPRTK)
  cepgen_build(CepGenExprtk
      SOURCES src/FunctionalExprTk.cpp
      EXT_HEADERS ${EXPRTK} # header-only
      INSTALL_COMPONENT exprtk)
endif()

#--- searching for TinyExpr
find_path(TINYEXPR_SRC NAMES tinyexpr.c HINTS $ENV{TINYEXPR_DIR})
find_path(TINYEXPR_INCLUDE NAMES tinyexpr.h HINTS $ENV{TINYEXPR_DIR})
if(TINYEXPR_SRC AND TINYEXPR_INCLUDE)
  cepgen_build(CepGenTinyExpr
      SOURCES src/FunctionalTinyExpr.cpp ${TINYEXPR_SRC}/tinyexpr.c
      EXT_HEADERS ${TINYEXPR_INCLUDE}
      INSTALL_COMPONENT tinyexpr)
endif()

#--- searching for libmatheval
find_library(LIBMATHEVAL NAMES matheval HINTS $ENV{LIBMATHEVAL_DIR} PATH_SUFFIXES lib)
find_path(LIBMATHEVAL_INCLUDE NAMES matheval.h HINTS $ENV{LIBMATHEVAL_DIR} PATH_SUFFIXES include)
if(LIBMATHEVAL)
  cepgen_build(CepGenLibmatheval
      SOURCES src/FunctionalLibmatheval.cpp
      EXT_LIBS ${LIBMATHEVAL}
      EXT_HEADERS ${LIBMATHEVAL_INCLUDE}
      INSTALL_COMPONENT libmatheval)
endif()
