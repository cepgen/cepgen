#--- searching for Cuba
find_path(CUBA_INCLUDE tools HINTS $ENV{CUBA_DIR})
find_library(CUBA cuba HINTS $ENV{CUBA_DIR})

if(NOT CUBA)
  return()
endif()

#----- build the object

cepgen_build(CepGenCuba SOURCES IntegratorCuba*.cpp
    EXT_LIBS ${CUBA}
    EXT_HEADERS ${CUBA_INCLUDE}
    INSTALL_COMPONENT cuba)

#set_target_properties(CepGenCore PROPERTIES COMPILE_FLAGS "-ffast-math -lm")
