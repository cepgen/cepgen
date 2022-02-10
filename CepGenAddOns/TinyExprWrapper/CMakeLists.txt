#--- searching for TinyExpr
find_path(TINYEXPR_SRC NAMES tinyexpr.c HINTS $ENV{TINYEXPR_DIR})
find_path(TINYEXPR_INCLUDE NAMES tinyexpr.h HINTS $ENV{TINYEXPR_DIR})

if(NOT TINYEXPR_SRC OR NOT TINYEXPR_INCLUDE)
  return()
endif()

file(GLOB sources *.cpp ${TINYEXPR_SRC}/tinyexpr.c)

#----- build the object

cepgen_build(CepGenTinyExpr SOURCES ${sources}
    EXT_HEADERS ${TINYEXPR_INCLUDE}
    INSTALL_COMPONENT tinyexpr)
