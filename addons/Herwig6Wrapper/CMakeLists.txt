#--- searching for Herwig 6
set(HERWIG6_DIRS ${CMAKE_CURRENT_BINARY_DIR} $ENV{HERWIG6_DIR} /usr /usr/local /opt/herwig6)
find_library(HERWIG6 NAMES herwig6 herwig herwig_pdfdummy HINTS ${HERWIG6_DIRS} PATH_SUFFIXES lib)

if(NOT HERWIG6)
  file(GLOB HERWIG6_SRC ${PROJECT_SOURCE_DIR}/External/herwig65*.f)
  file(GLOB HERWIG6_INC ${PROJECT_SOURCE_DIR}/External)
  if(EXISTS $ENV{HERWIG6_SRC})
    file(GLOB herwig6_src $ENV{HERWIG6_SRC})
  elseif(EXISTS ${HERWIG6_SRC})
    file(GLOB herwig6_src ${HERWIG6_SRC})
  else()
    return()
  endif()
  message(STATUS "Herwig 6 source found in ${herwig6_src}")
  add_library(herwig6 SHARED ${herwig6_src})
  target_include_directories(herwig6 PUBLIC ${HERWIG6_INC})
  install(TARGETS herwig6 DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT herwig6)
  set(HERWIG6 herwig6)
endif()

#----- build the object
enable_testing()
cepgen_build(CepGenHerwig6
    SOURCES src/*.cpp
    LIBRARIES ${HERWIG6}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    OPTIONS "-Wno-misleading-indentation"
    COMPONENT herwig6)
install(DIRECTORY CepGenHerwig6
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT herwig6)
