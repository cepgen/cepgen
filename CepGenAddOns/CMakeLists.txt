file(GLOB io_sources OutputModules/TextHandler.cpp)
set(hadr_sources)

include_directories(${PROJECT_SOURCE_DIR})

set(addons_libraries)
set(addons_definitions)
set(addons_includes)

#----- check the external dependencies for physics utilities

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
  list(APPEND io_sources EventInterfaces/PythiaEventInterface.cpp OutputModules/LHEFPythiaHandler.cpp)
endif()

#--- linking with ROOT/Delphes

if(ROOT_FOUND)
  message(STATUS "ROOT found in ${ROOT_LIBRARY_DIR}")
  list(APPEND addons_libraries ${ROOT_LIBRARIES})
  list(APPEND addons_definitions ${ROOT_DEFINITIONS})
  list(APPEND addons_includes ${ROOT_INCLUDE_DIRS})
  list(APPEND io_sources EventInterfaces/ROOTTreeInfo.cpp OutputModules/ROOTHistsHandler.cpp OutputModules/ROOTTreeHandler.cpp)
  if(DELPHES)
    message(STATUS "Delphes found in ${DELPHES}")
    list(APPEND addons_libraries ${DELPHES})
    list(APPEND addons_includes ${DELPHES_INCLUDE} ${DELPHES_EXT_INCLUDE})
    list(APPEND io_sources OutputModules/DelphesHandler.cpp)
  endif()
endif()

#--- linking with HepMC

if(HEPMC_LIB)
  message(STATUS "HepMC found in ${HEPMC_LIB}")
  message(STATUS "HepMC headers found in ${HEPMC_INCLUDE}")
  if(HEPMC_LIB MATCHES ".*HepMC3.?.so")
    message(STATUS "HepMC version 3 found")
    if(HEPMC_ROOT_LIB)
      message(STATUS "HepMC ROOT I/O library found")
      list(APPEND addons_libraries ${HEPMC_ROOT_LIB})
      list(APPEND addons_definitions -DHEPMC3_ROOTIO)
    endif()
    if(HEPMC_EXTRA_PLUGINS)
      message(STATUS "HepMC example plugins found in ${HEPMC_EXTRA_PLUGINS}")
      file(GLOB ConvertExample_sources ${HEPMC_EXTRA_PLUGINS}/ConvertExample/src/*.cc)
      set(hepmc_example_includes ${HEPMC_EXTRA_PLUGINS}/ConvertExample/include)
      if(ROOT_FOUND)
        list(APPEND hepmc_example_includes ${ROOT_INCLUDE_DIRS})
      else()
        list(REMOVE_ITEM ConvertExample_sources ${HEPMC_EXTRA_PLUGINS}/ConvertExample/src/WriterRootTreeOPAL.cc)
      endif()
      add_library(HepMCExamples SHARED ${ConvertExample_sources})
      target_include_directories(HepMCExamples PUBLIC ${hepmc_example_includes})
      if(ROOT_FOUND)
        set_property(TARGET HepMCExamples PROPERTY CXX_STANDARD 14)
      endif()
      list(APPEND addons_libraries HepMCExamples)
      list(APPEND addons_includes ${HEPMC_EXTRA_PLUGINS})
      list(APPEND addons_definitions -DHEPMC3_EXTRA_PLUGINS)
    endif()
    list(APPEND addons_definitions -DHEPMC3)
  endif()
  list(APPEND io_sources EventInterfaces/HepMCEventInterface.cpp OutputModules/HepMCHandler.cpp)
  if(NOT PYTHIA8)
    message(STATUS "HepMC will be used for LHEF output")
    list(APPEND io_sources OutputModules/LHEFHepMCHandler.cpp)
  endif()
  list(APPEND addons_libraries ${HEPMC_LIB})
  list(APPEND addons_includes ${HEPMC_INCLUDE})
endif()

#--- linking with ProMC

if(PROMC_LIB)
  message(STATUS "ProMC found in ${PROMC_LIB}")
  list(APPEND addons_libraries ${PROMC_LIB})
  list(APPEND addons_includes ${PROMC_INCLUDE} ${PROMC_EXT_INCLUDE})
  list(APPEND io_sources OutputModules/ProMCHandler.cpp)
endif()

#--- linking with YODA

if(YODA)
  message(STATUS "YODA found in ${YODA}")
  list(APPEND addons_libraries ${YODA})
  list(APPEND addons_includes ${YODA_INCLUDE})
  list(APPEND io_sources OutputModules/YODAHistsHandler.cpp)
endif()

#----- build the objects

add_library(CepGenAddOns SHARED ${io_sources} ${hadr_sources})
target_link_libraries(CepGenAddOns ${addons_libraries} CepGenEvent)
target_compile_definitions(CepGenAddOns PUBLIC ${addons_definitions})
target_include_directories(CepGenAddOns PUBLIC ${addons_includes})
if(ROOT_FOUND)
  set_property(TARGET CepGenAddOns PROPERTY CXX_STANDARD 14)
endif()

#----- installation rules

install(TARGETS CepGenAddOns DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)

