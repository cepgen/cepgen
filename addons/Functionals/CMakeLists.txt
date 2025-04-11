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

#--- searching for fparser
find_path(FunctionParser_INCLUDE NAMES fparser.hh HINTS $ENV{FunctionParser_DIR} PATH_SUFFIXES include)
find_path(FunctionParser_SOURCES NAMES fparser.cc HINTS $ENV{FunctionParser_DIR} PATH_SUFFIXES include)
if(FunctionParser_INCLUDE AND FunctionParser_SOURCES)
  file(GLOB fparser_SRC ${FunctionParser_SOURCES}/*.cc)
  cepgen_build(CepGenFunctionParser
      SOURCES src/FunctionalFunctionParser.cpp ${fparser_SRC}
      INCLUDES ${FunctionParser_INCLUDE}
      COMPONENT fparser)
endif()
