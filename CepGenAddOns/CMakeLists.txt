include_directories(${PROJECT_SOURCE_DIR})

file(GLOB output_sources OutputModules/EventDump.cpp OutputModules/TextHandler.cpp)
set(modif_sources)
set(addons_libraries)
set(addons_definitions)
set(addons_includes)

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

#----- build the objects

add_library(CepGenAddOns SHARED ${modif_sources} ${output_sources})
target_compile_definitions(CepGenAddOns PRIVATE ${addons_definitions})
target_include_directories(CepGenAddOns PRIVATE ${addons_includes})
target_link_libraries(CepGenAddOns PRIVATE ${addons_libraries})

#----- propagate global variables to parent scope

set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)
set(CEPGEN_HAS_FUNCTIONALS ${CEPGEN_HAS_FUNCTIONALS} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenAddOns DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lib)
