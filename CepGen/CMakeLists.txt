file(GLOB core_sources Core/*.cpp *.cpp)
file(GLOB phys_sources Physics/*.cpp)
file(GLOB sf_sources StructureFunctions/*.cpp)
file(GLOB io_sources IO/GenericExportHandler.cpp)
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

#--- linking with Pythia 8

if(PYTHIA8)
  message(STATUS "Pythia 8 found in ${PYTHIA8}")
  message(STATUS "Pythia 8 will be used for LHEF output")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
  include_directories(${PYTHIA8_INCLUDE})
  add_definitions(-DPYTHIA8)
  list(APPEND hadr_sources "Hadronisers/Pythia8Hadroniser.cpp")
  list(APPEND hadr_sources "Hadronisers/PythiaEventInterface.cpp")
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

#--- linking with Delphes

if(DELPHES)
  message(STATUS "Delphes found in ${DELPHES}")
  find_package(ROOT COMPONENTS Core QUIET)
  if(ROOT_FOUND)
    message(STATUS "ROOT found in ${ROOT_INCLUDE_DIRS}")
    list(APPEND addons_libraries ${DELPHES_LIB} ${ROOT_LIBRARIES} lzma)
    list(APPEND io_sources "IO/DelphesHandler.cpp")
    include_directories(${DELPHES_INCLUDE})
    include_directories(${ROOT_INCLUDE_DIRS})
  else()
    message(STATUS "ROOT not found! Disabling Delphes module.")
  endif()
endif()

#----- build the objects

add_library(CepGenCore SHARED ${core_sources} ${phys_sources} ${sf_sources})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_CORE_REQS})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_STRF_REQS})

add_library(CepGenAddOns SHARED ${io_sources} ${hadr_sources})
target_link_libraries(CepGenAddOns ${addons_libraries})
target_link_libraries(CepGenAddOns CepGenEvent)

#----- installation rules

install(TARGETS CepGenCore DESTINATION lib)
install(TARGETS CepGenAddOns DESTINATION lib)

