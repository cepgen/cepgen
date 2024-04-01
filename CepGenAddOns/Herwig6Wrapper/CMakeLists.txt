#--- searching for Herwig 6
set(HERWIG6_DIRS $ENV{HERWIG6_DIR} /usr /usr/local /opt/herwig6)
find_library(HERWIG6 herwig6 HINTS ${HERWIG6_DIRS} PATH_SUFFIXES lib)

if(NOT HERWIG6)
  file(GLOB HERWIG6_SRC ${PROJECT_SOURCE_DIR}/External/herwig65*.f)
  if(EXISTS $ENV{HERWIG6_SRC})
    file(GLOB herwig6_src $ENV{HERWIG6_SRC})
  elseif(EXISTS ${HERWIG6_SRC})
    file(GLOB herwig6_src ${HERWIG6_SRC})
  else()
    return()
  endif()
  message(STATUS "Herwig 6 source found in ${herwig6_src}")
  add_library(herwig6 SHARED ${herwig6_src})
  install(TARGETS herwig6 DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT herwig6)
  target_compile_options(herwig6 PRIVATE "-Wno-tabs -Wno-maybe-uninitialized -Wno-integer-division -Wno-unused-variable -Wno-unused-dummy-argument")
  set(HERWIG6 herwig6)
endif()

#----- build the objects

enable_testing()
cepgen_build(CepGenHerwig6 SOURCES *.cpp
    EXT_LIBS ${HERWIG6}
    OPTIONS "-Wno-misleading-indentation"
    INSTALL_COMPONENT herwig6)
