cepgen_build(CepGenAddOns
    SOURCES OutputModules/EventDump.cpp OutputModules/TextHandler.cpp
    INSTALL_COMPONENT addons)

#----- check the external dependencies for physics utilities

add_subdirectory(BoostWrapper)
add_subdirectory(Pythia6Wrapper)
add_subdirectory(Pythia8Wrapper)
add_subdirectory(ROOTWrapper)
add_subdirectory(HepMCWrapper)
add_subdirectory(ProMCWrapper)
add_subdirectory(LHAPDFWrapper)
add_subdirectory(RivetWrapper)
add_subdirectory(APFELWrapper)
add_subdirectory(MadGraphWrapper)

#----- propagate global variables to parent scope

set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)
set(CEPGEN_HAS_FUNCTIONALS ${CEPGEN_HAS_FUNCTIONALS} PARENT_SCOPE)

