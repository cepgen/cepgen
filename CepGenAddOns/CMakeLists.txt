#----- check the external dependencies for physics utilities

if(CMAKE_BUILD_TESTS)
  enable_testing()
endif()

add_subdirectory(APFELWrapper)
add_subdirectory(BoostWrapper)
add_subdirectory(CubaWrapper)
add_subdirectory(DelphesWrapper)
add_subdirectory(HepMC2Wrapper)
add_subdirectory(HepMC3Wrapper)
add_subdirectory(LHAPDFWrapper)
add_subdirectory(MadGraphWrapper)
add_subdirectory(MatplotlibWrapper)
add_subdirectory(ProMCWrapper)
add_subdirectory(Pythia6Wrapper)
add_subdirectory(Pythia8Wrapper)
add_subdirectory(PythonWrapper)
add_subdirectory(RivetWrapper)
add_subdirectory(ROOTWrapper)
add_subdirectory(TauolaWrapper)
add_subdirectory(TopdrawerWrapper)

add_subdirectory(Functionals)

#----- propagate global variables to parent scope

set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

