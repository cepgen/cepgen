if(HEPMC_LIB)
  message(STATUS "HepMC found in ${HEPMC_LIB}; headers found in ${HEPMC_INCLUDE}")
else()
  message(FATAL_ERROR "Trying to build HepMC wrapper library while no HepMC installation was found on system!")
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources *.cpp)

#----- build the objects

add_library(CepGenHepMC SHARED ${sources})
if(HEPMC_LIB MATCHES ".*HepMC3.?.so")
  message(STATUS "HepMC version 3 found")
  set(compile_def -DHEPMC3)
  if(HEPMC_ROOT_LIB)
    message(STATUS "HepMC ROOT I/O library found")
    list(APPEND addons_libraries ${HEPMC_ROOT_LIB})
    list(APPEND compile_def -DHEPMC3_ROOTIO)
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
    list(APPEND HEPMC_LIB HepMCExamples)
    list(APPEND HEPMC_INCLUDE ${HEPMC_EXTRA_PLUGINS})
    list(APPEND compile_def -DHEPMC3_EXTRA_PLUGINS)
  endif()
  target_compile_definitions(CepGenHepMC PUBLIC ${compile_def})
endif()
target_link_libraries(CepGenHepMC PUBLIC ${HEPMC_LIB})
target_include_directories(CepGenHepMC PUBLIC ${HEPMC_INCLUDE})
list(APPEND CEPGEN_ADDONS_EXT CepGenHepMC)
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenHepMC DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT hepmc)


