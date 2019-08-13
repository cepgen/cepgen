file(GLOB core_sources Core/*.cpp *.cpp)
file(GLOB phys_sources Physics/*.cpp)
file(GLOB sf_sources StructureFunctions/*.cpp)
file(GLOB io_sources IO/GenericExportHandler.cpp IO/TextHandler.cpp)
file(GLOB hadr_sources Hadronisers/GenericHadroniser.cpp)

include_directories(${PROJECT_SOURCE_DIR})

message(STATUS "GSL found in ${GSL_LIB}")

set(core_libraries ${GSL_LIB} ${GSL_CBLAS_LIB})
set(core_definitions)
set(core_includes)
set(addons_libraries)
set(addons_definitions)
set(addons_includes)

#----- check the external dependencies for SFs

if(GRV_SRC)
  message(STATUS "GRV PDFset found in ${GRV_SRC}!")
  list(APPEND sf_sources ${GRV_SRC})
  list(APPEND core_definitions -DGRVPDF)
endif()

#----- check the external dependencies for physics utilities

if(ALPHAS_SRC)
  message(STATUS "alphaS evolution algorithm found in ${ALPHAS_SRC}")
  list(APPEND phys_sources ${ALPHAS_SRC})
  list(APPEND core_definitions -DALPHA_S)
endif()
if(MUPARSER)
  message(STATUS "Formulas parsing will use muParser found in ${MUPARSER}")
  list(APPEND core_libraries ${MUPARSER})
elseif(EXPRTK)
  message(STATUS "Formulas parsing will use exprtk found in ${EXPRTK}")
elseif(ROOT_FOUND)
  message(STATUS "Formulas parsing will use ROOT")
  list(APPEND core_libraries ${ROOT_LIBRARIES})
endif()

#--- linking with LHAPDF

if(LHAPDF)
  message(STATUS "LHAPDF found in ${LHAPDF}")
  list(APPEND core_libraries ${LHAPDF})
  list(APPEND core_includes ${LHAPDF_INCLUDE})
else()
  file(GLOB partonic_sf StructureFunctions/Partonic.cpp)
  list(REMOVE_ITEM sf_sources ${partonic_sf})
endif()

#--- linking with Pythia 6

if(PYTHIA6)
  message(STATUS "Pythia 6 found in ${PYTHIA6}")
  list(APPEND hadr_sources Hadronisers/Pythia6Hadroniser.cpp)
  list(APPEND addons_libraries ${PYTHIA6})
  if(PYTHIA6DUMMY)
    message(STATUS "Pythia 6 addons found in ${PYTHIA6DUMMY}")
    list(APPEND addons_libraries ${PYTHIA6DUMMY})
  endif()
elseif(EXISTS $ENV{PYTHIA6_SRC})
  file(GLOB pythia6_src $ENV{PYTHIA6_SRC})
  message(STATUS "Pythia 6 source found in ${pythia6_src}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-tabs -Wno-maybe-uninitialized -Wno-integer-division -Wno-unused-variable -Wno-unused-dummy-argument")
  add_library(pythia6 SHARED ${pythia6_src})
  list(APPEND hadr_sources Hadronisers/Pythia6Hadroniser.cpp)
  list(APPEND addons_libraries pythia6)
endif()

#--- linking with Pythia 8

if(PYTHIA8)
  message(STATUS "Pythia 8 found in ${PYTHIA8}")
  message(STATUS "Pythia 8 will be used for LHEF output")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
  list(APPEND addons_libraries ${PYTHIA8} dl)
  list(APPEND addons_includes ${PYTHIA8_INCLUDE})
  list(APPEND addons_definitions -DPYTHIA8)
  list(APPEND hadr_sources Hadronisers/Pythia8Hadroniser.cpp)
  list(APPEND io_sources IO/PythiaEventInterface.cpp IO/LHEFPythiaHandler.cpp)
endif()

#--- linking with HepMC

if(HEPMC_LIB)
  message(STATUS "HepMC found in ${HEPMC_LIB}")
  if(HEPMC_LIB MATCHES ".*HepMC3.?.so")
    message(STATUS "HepMC version 3 found")
    if(HEPMC_ROOT_LIB)
      message(STATUS "HepMC ROOT I/O library found")
      list(APPEND addons_libraries ${HEPMC_ROOT_LIB})
      list(APPEND addons_definitions -DHEPMC3_ROOTIO)
    endif()
    list(APPEND addons_definitions -DHEPMC3)
  endif()
  list(APPEND io_sources IO/HepMCHandler.cpp IO/HepMCEventInterface.cpp)
  if(NOT PYTHIA8)
    message(STATUS "HepMC will be used for LHEF output")
    list(APPEND io_sources IO/LHEFHepMCHandler.cpp)
  endif()
  list(APPEND addons_libraries ${HEPMC_LIB})
  list(APPEND addons_includes ${HEPMC_INCLUDE})
endif()

#--- linking with ProMC

if(PROMC_LIB)
  message(STATUS "ProMC found in ${PROMC_LIB}")
  list(APPEND addons_libraries ${PROMC_LIB})
  list(APPEND addons_includes ${PROMC_INCLUDE} ${PROMC_EXT_INCLUDE})
  list(APPEND io_sources IO/ProMCHandler.cpp)
endif()

#--- linking with ROOT/Delphes

if(ROOT_FOUND)
  message(STATUS "ROOT found in ${ROOT_LIBRARY_DIR}")
  list(APPEND addons_libraries ${ROOT_LIBRARIES})
  list(APPEND addons_definitions ${ROOT_DEFINITIONS})
  list(APPEND addons_includes ${ROOT_INCLUDE_DIRS})
  list(APPEND io_sources IO/ROOTHistsHandler.cpp IO/ROOTTreeHandler.cpp)
  if(DELPHES)
    message(STATUS "Delphes found in ${DELPHES}")
    list(APPEND addons_libraries ${DELPHES})
    list(APPEND addons_includes ${DELPHES_INCLUDE} ${DELPHES_EXT_INCLUDE})
    list(APPEND io_sources IO/DelphesHandler.cpp)
  endif()
endif()

#--- linking with YODA

if(YODA)
  message(STATUS "YODA found in ${YODA}")
  list(APPEND addons_libraries ${YODA})
  list(APPEND addons_includes ${YODA_INCLUDE})
  list(APPEND io_sources IO/YODAHistsHandler.cpp)
endif()

#----- build the objects

add_library(CepGenCore SHARED ${core_sources} ${phys_sources} ${sf_sources})
target_link_libraries(CepGenCore ${core_libraries})
target_compile_definitions(CepGenCore PUBLIC ${core_definitions})
target_include_directories(CepGenCore PUBLIC ${core_includes})

add_library(CepGenAddOns SHARED ${io_sources} ${hadr_sources})
target_link_libraries(CepGenAddOns ${addons_libraries} CepGenEvent)
target_compile_definitions(CepGenAddOns PUBLIC ${addons_definitions})
target_include_directories(CepGenAddOns PUBLIC ${addons_includes})

#----- installation rules

install(TARGETS CepGenCore DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)
install(TARGETS CepGenAddOns DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)

