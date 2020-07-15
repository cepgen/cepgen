include_directories(${PROJECT_SOURCE_DIR})

file(GLOB output_sources OutputModules/EventDump.cpp OutputModules/TextHandler.cpp)
set(modif_sources)
set(addons_libraries)
set(addons_definitions)
set(addons_includes)

#----- check the external dependencies for physics utilities

if(Boost_FOUND)
  add_subdirectory(BoostWrapper)
endif()
add_subdirectory(Pythia6Wrapper)
if(PYTHIA8)
  add_subdirectory(Pythia8Wrapper)
endif()
if(ROOT_FOUND)
  add_subdirectory(ROOTWrapper)
endif()
if(HEPMC_LIB)
  add_subdirectory(HepMCWrapper)
endif()
if(LHAPDF)
  add_subdirectory(LHAPDFWrapper)
endif()
if(YODA)
  add_subdirectory(RivetWrapper)
endif()

#--- linking with ROOT/Delphes

if(DELPHES)
  message(STATUS "Delphes found in ${DELPHES}")
  list(APPEND addons_libraries ${DELPHES})
  list(APPEND addons_includes ${DELPHES_INCLUDE} ${DELPHES_EXT_INCLUDE})
  list(APPEND output_sources OutputModules/DelphesHandler.cpp)
endif()

#--- linking with ProMC

if(PROMC_LIB)
  message(STATUS "ProMC found in ${PROMC_LIB}")
  list(APPEND addons_libraries ${PROMC_LIB})
  list(APPEND addons_includes ${PROMC_INCLUDE} ${PROMC_EXT_INCLUDE})
  list(APPEND output_sources OutputModules/ProMCHandler.cpp)
endif()

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
