set(SHERPA_DIRS $ENV{SHERPA_DIR} $ENV{SHERPA_ROOT} /usr /usr/local /opt/sherpa)
find_package(SHERPA-MC 3.0 QUIET CONFIG HINTS ${SHERPA_DIRS} PATH_SUFFIXES share/SHERPA-MC lib/cmake/SHERPA-MC)
if(NOT SHERPA-MC_FOUND)
  return()
endif()
message(STATUS "Sherpa ${SHERPA-MC_VERSION} found in ${SHERPA-MC_DIR}")

find_path(SHERPA_INCLUDE SHERPA/Main/Sherpa.H
    HINTS ${SHERPA-MC_DIR}/../.. ${SHERPA-MC_DIR}/../../.. ${SHERPA_DIRS}
    PATH_SUFFIXES include/SHERPA-MC include)
if(NOT SHERPA_INCLUDE)
  return()
endif()

#----- build the object
cepgen_build(CepGenSherpa
    SOURCES src/*.cpp
    LIBRARIES SHERPA::All
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${SHERPA_INCLUDE}
    COMPONENT sherpa)
if(NOT APPLE)
  target_link_options(CepGenSherpa PRIVATE "LINKER:--no-as-needed")
endif()
cpack_add_component(sherpa
    DISPLAY_NAME "CepGen Sherpa wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Sherpa fragmentation and hadron decays algorithms"
    DEPENDS lib)
install(DIRECTORY CepGenSherpa
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT sherpa)
