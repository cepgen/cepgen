set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_SET_DESTDIR TRUE)
set(CPACK_PACKAGE_NAME "CepGen")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A generic central exclusive processes event generator")
set(CPACK_PACKAGE_RELOCATABLE FALSE)
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "Laurent Forthomme <cepgen@projects.hepforge.org>")
set(CPACK_PACKAGE_LICENSE GPL3+)
set(CPACK_PACKAGE_VENDOR "Helsinki Institute of Physics")
set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/doc/small-cepgen-logo.png")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CepGen/README")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
find_program(RPMBUILD rpmbuild QUIET)
if(RPMBUILD)
    #--- RPM information
    message(STATUS "RPM packaging utilitaries detected")
    set(CPACK_GENERATOR "RPM")
    set(CPACK_RPM_MAIN_COMPONENT lib)
    set(CPACK_RPM_FILE_NAME RPM-DEFAULT)
    set(CPACK_RPM_CHANGELOG_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CHANGELOG")
    set(CPACK_RPM_PACKAGE_URL https://cepgen.hepforge.org/)
    #set(CPACK_RPM_PACKAGE_SOURCE TRUE)
    set(CPACK_RPM_ENABLE_COMPONENT_DEPENDS TRUE)
    #--- packages interdependencies
    set(CPACK_RPM_LIB_PACKAGE_REQUIRES "gsl >= 2.1, python >= 2.7")
    set(CEPGEN_MIN_REQ "cepgen == ${VERSION}")
    set(CPACK_RPM_DEVEL_PACKAGE_REQUIRES ${CEPGEN_MIN_REQ})
    set(CPACK_RPM_DEVEL_PACKAGE_ARCHITECTURE noarch)
    set(CPACK_RPM_ROOT_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, root >= 6.0")
    set(CPACK_RPM_PYTHIA8_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, pythia8 >= 8.2.30")
    set(CPACK_RPM_BOOST_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, boost >= 1.33")
    set(CPACK_RPM_LHAPDF_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, lhapdf")
    set(CPACK_RPM_HEPMC_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, HepMC >= 2.01")
    set(CPACK_RPM_RIVET_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}")
    #set(CPACK_RPM_PACKAGE_AUTOREQPROV " no")
    list(APPEND CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST /usr)
    list(APPEND CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST /usr/bin)
    list(APPEND CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST /usr/lib)
    list(APPEND CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST /usr/lib64)
    list(APPEND CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST /usr/include)
else()
    #--- DEB information
    message(STATUS "DEB packaging may be performed")
    set(CPACK_GENERATOR "DEB")
    set(CPACK_DEB_MAIN_COMPONENT lib)
    set(CPACK_DEB_CHANGELOG_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CHANGELOG")
    set(CPACK_DEB_PACKAGE_URL https://cepgen.hepforge.org/)
    set(CPACK_DEB_ENABLE_COMPONENT_DEPENDS TRUE)
    #--- packages interdependencies
    set(CPACK_DEB_LIB_PACKAGE_REQUIRES "gsl >= 2.1, python >= 2.7")
    set(CEPGEN_MIN_REQ "cepgen == ${VERSION}")
    set(CPACK_DEB_DEVEL_PACKAGE_REQUIRES ${CEPGEN_MIN_REQ})
    set(CPACK_DEB_DEVEL_PACKAGE_ARCHITECTURE noarch)
    set(CPACK_DEB_ROOT_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, root >= 6.0")
    set(CPACK_DEB_PYTHIA8_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, pythia8 >= 8.2.30")
    set(CPACK_DEB_BOOST_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, boost >= 1.33")
    set(CPACK_DEB_LHAPDF_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, lhapdf")
    set(CPACK_DEB_HEPMC_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}, HepMC >= 2.01")
    set(CPACK_DEB_RIVET_PACKAGE_REQUIRES "${CEPGEN_MIN_REQ}")
endif()
include(CPack)
#--- register packages
cpack_add_component(lib
    DISPLAY_NAME "CepGen core library"
    DESCRIPTION "The full set of core libraries embedded within CepGen"
    REQUIRED)
#set(CPACK_COMPONENTS_ALL lib devel root pythia8 boost lhapdf)
cpack_add_component(devel
    DISPLAY_NAME "CepGen development headers"
    DESCRIPTION "Collection of C and C++ includes for the development of CepGen-dependent libraries"
    DEPENDS lib)
cpack_add_component(root
    DISPLAY_NAME "CepGen ROOT wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the ROOT library"
    DEPENDS lib)
cpack_add_component(pythia8
    DISPLAY_NAME "CepGen Pythia 8 wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to Pythia 8"
    DEPENDS lib)
cpack_add_component(boost
    DISPLAY_NAME "CepGen Boost wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Boost library"
    DEPENDS lib)
cpack_add_component(hepmc
    DISPLAY_NAME "CepGen HepMC wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the HepMC library"
    DEPENDS lib)
cpack_add_component(lhapdf
    DISPLAY_NAME "CepGen LHAPDF wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the LHAPDF library"
    DEPENDS lib)
cpack_add_component(rivet
    DISPLAY_NAME "CepGen Rivet wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Rivet library"
    DEPENDS lib)
#message(STATUS ">>> ${CPACK_COMPONENTS_ALL}")
