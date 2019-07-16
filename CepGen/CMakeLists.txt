file(GLOB core_sources Core/*.cpp *.cpp)
file(GLOB phys_sources Physics/*.cpp)
file(GLOB sf_sources StructureFunctions/*.cpp)
file(GLOB io_sources IO/GenericExportHandler.cpp IO/TextHandler.cpp)
file(GLOB hadr_sources Hadronisers/GenericHadroniser.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- check the external dependencies for SFs

set(GRV_PATH ${PROJECT_SOURCE_DIR}/External)
file(GLOB grv_sources ${GRV_PATH}/grv_*.f)
if(grv_sources)
  message(STATUS "GRV PDFset found in ${grv_sources}!")
  add_definitions(-DGRVPDF)
  list(APPEND sf_sources ${grv_sources})
else()
  message(STATUS "GRV PDFset not found. Will proceed without it")
endif()

#----- check the external dependencies for physics utilities

if(alphas_sources)
  list(APPEND phys_sources ${alphas_sources})
endif()

set(addons_libraries "")
set(strf_libraries ${GSL_LIB} ${GSL_CBLAS_LIB})

#--- linking with LHAPDF

if(LHAPDF)
  message(STATUS "LHAPDF found in ${LHAPDF}")
  list(APPEND strf_libraries ${LHAPDF})
  include_directories(${LHAPDF_INCLUDE})
else()
  file(GLOB partonic_sf StructureFunctions/Partonic.cpp)
  list(REMOVE_ITEM sf_sources ${partonic_sf})
endif()

#--- linking with Pythia 6

if(PYTHIA6)
  message(STATUS "Pythia 6 found in ${PYTHIA6}")
  list(APPEND hadr_sources "Hadronisers/Pythia6Hadroniser.cpp")
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
  list(APPEND hadr_sources "Hadronisers/Pythia6Hadroniser.cpp")
  list(APPEND addons_libraries pythia6)
endif()

#--- linking with Pythia 8

if(PYTHIA8)
  message(STATUS "Pythia 8 found in ${PYTHIA8}")
  message(STATUS "Pythia 8 will be used for LHEF output")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
  include_directories(${PYTHIA8_INCLUDE})
  add_definitions(-DPYTHIA8)
  list(APPEND hadr_sources "Hadronisers/Pythia8Hadroniser.cpp")
  list(APPEND io_sources "IO/PythiaEventInterface.cpp")
  list(APPEND io_sources "IO/LHEFPythiaHandler.cpp")
  list(APPEND addons_libraries ${PYTHIA8} dl)
endif()

#--- linking with HepMC

if(HEPMC_LIB)
  message(STATUS "HepMC found in ${HEPMC_LIB}")
  if(HEPMC_LIB MATCHES ".*HepMC3.?.so")
    message(STATUS "HepMC version 3 found")
    if(HEPMC_ROOT_LIB)
      message(STATUS "HepMC ROOT I/O library found")
      add_definitions(-DHEPMC3_ROOTIO)
      list(APPEND addons_libraries ${HEPMC_ROOT_LIB})
    endif()
    add_definitions(-DHEPMC3)
  endif()
  list(APPEND io_sources "IO/HepMCHandler.cpp")
  if(NOT PYTHIA8)
    message(STATUS "HepMC will be used for LHEF output")
    list(APPEND io_sources "IO/LHEFHepMCHandler.cpp")
  endif()
  list(APPEND addons_libraries ${HEPMC_LIB})
  include_directories(${HEPMC_INCLUDE})
endif()

#--- linking with ROOT/Delphes

if(ROOT_FOUND)
  message(STATUS "ROOT found in ${ROOT_LIBRARY_DIR}")
  list(APPEND addons_libraries ${ROOT_LIBRARIES})
  list(APPEND io_sources "IO/ROOTTreeHandler.cpp")
  include_directories(${ROOT_INCLUDE_DIRS})
  link_directories(${ROOT_LIBRARY_DIR})
  if(DELPHES)
    message(STATUS "Delphes found in ${DELPHES}")
    list(APPEND addons_libraries ${DELPHES} ${TBB})
    list(APPEND io_sources "IO/DelphesHandler.cpp")
    include_directories(${DELPHES_INCLUDE})
    message(STATUS ${DELPHES_INCLUDE})
  endif()
endif()

#----- build the objects

add_library(CepGenCore SHARED ${core_sources} ${phys_sources} ${sf_sources})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_CORE_REQS})
target_link_libraries(CepGenCore ${strf_libraries})

add_library(CepGenAddOns SHARED ${io_sources} ${hadr_sources})
target_link_libraries(CepGenAddOns ${addons_libraries})
target_link_libraries(CepGenAddOns CepGenEvent)

#----- installation rules

install(TARGETS CepGenCore DESTINATION lib)
install(TARGETS CepGenAddOns DESTINATION lib)

