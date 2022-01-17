#----- check the external dependencies for physics utilities

add_subdirectory(APFELWrapper)
add_subdirectory(BoostWrapper)
add_subdirectory(DelphesWrapper)
add_subdirectory(HepMC2Wrapper)
add_subdirectory(HepMC3Wrapper)
add_subdirectory(LHAPDFWrapper)
add_subdirectory(MadGraphWrapper)
add_subdirectory(ProMCWrapper)
add_subdirectory(Pythia6Wrapper)
add_subdirectory(Pythia8Wrapper)
add_subdirectory(RivetWrapper)
add_subdirectory(ROOTWrapper)

add_subdirectory(Functionals)

#----- propagate global variables to parent scope

set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

