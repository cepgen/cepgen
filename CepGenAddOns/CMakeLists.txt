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
if(PYTHIA8)
  add_subdirectory(PythiaWrapper)
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

#--- linking with Pythia 6

if(PYTHIA6)
  message(STATUS "Pythia 6 found in ${PYTHIA6}")
  list(APPEND modif_sources Hadronisers/Pythia6Hadroniser.cpp)
  list(APPEND addons_libraries ${PYTHIA6})
  if(PYTHIA6DUMMY)
    message(STATUS "Pythia 6 addons found in ${PYTHIA6DUMMY}")
    list(APPEND addons_libraries ${PYTHIA6DUMMY})
  endif()
elseif(EXISTS $ENV{PYTHIA6_SRC})
  file(GLOB pythia6_src $ENV{PYTHIA6_SRC})
  message(STATUS "Pythia 6 source found in ${pythia6_src}")
  add_library(pythia6 SHARED ${pythia6_src})
  set_property(TARGET pythia6 PROPERTY Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-tabs -Wno-maybe-uninitialized -Wno-integer-division -Wno-unused-variable -Wno-unused-dummy-argument")
  list(APPEND modif_sources Hadronisers/Pythia6Hadroniser.cpp)
  list(APPEND addons_libraries pythia6)
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
target_compile_definitions(CepGenAddOns PUBLIC ${addons_definitions})
target_include_directories(CepGenAddOns PUBLIC ${addons_includes})
target_link_libraries(CepGenAddOns PUBLIC ${addons_libraries})
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)
set(CEPGEN_HAS_FUNCTIONALS ${CEPGEN_HAS_FUNCTIONALS} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenAddOns DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lib)
