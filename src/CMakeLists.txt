file(GLOB core_sources
        Cards/*.cpp
        Core/*.cpp
        EventFilter/*.cpp
        Integration/*.cpp
        Modules/*.cpp
        OutputModules/*.cpp
        Utils/*.cpp
        *.cpp)
file(GLOB phys_sources
        Event/*.cpp
        FormFactors/*.cpp
        PartonFluxes/*.cpp
        Physics/*.cpp
        StructureFunctions/*.cpp)
file(GLOB proc_sources
        Process/*.cpp
        Process/*.f)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

#----- check the external dependencies for SFs/alphaS

file(GLOB GRV_SRC ${PROJECT_SOURCE_DIR}/External/grv_*.f)
if(GRV_SRC)
  message(STATUS "GRV PDFset found in ${GRV_SRC}!")
  list(APPEND phys_sources ${GRV_SRC})
else()
  message(STATUS "GRV PDFset not found. Will not build Szczurek-Uleshchenko SFs")
  file(GLOB su_sf StructureFunctions/SzczurekUleshchenko.cpp)
  list(REMOVE_ITEM phys_sources ${su_sf})
endif()

file(GLOB ALPHAS_SRC ${PROJECT_SOURCE_DIR}/External/alphaS.f)
if(ALPHAS_SRC)
  message(STATUS "alphaS evolution algorithm found in ${ALPHAS_SRC}")
  list(APPEND phys_sources ${ALPHAS_SRC})
else()
  file(GLOB pegasus_alphas Physics/AlphaSPEGASUS.cpp)
  file(REMOVE_ITEM phys_sources ${pegasus_alphas})
endif()

#----- build the objects

add_version_definition(${core_sources})

cepgen_build(CepGen
        SOURCES ${core_sources} ${phys_sources} ${proc_sources}
        LIBRARIES ${CEPGEN_CORE_EXT}
        COMPONENT lib)

file(GLOB exec_sources *.cc)
message(STATUS "Executables will be linked against ${CEPGEN_LIBRARIES}")
foreach(exec_src ${exec_sources})
  get_filename_component(exec_bin ${exec_src} NAME_WE)
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES})
  install(TARGETS ${exec_bin} DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT lib)
endforeach()

add_test(NAME cepgen_list COMMAND cepgen -l)
set_property(TEST cepgen_list PROPERTY ENVIRONMENT "CEPGEN_PATH=${CEPGEN_PATH}")
add_test(NAME cepgen_lpair COMMAND cepgen -i Cards/lpair_cfg.py -n 5000)
set_property(TEST cepgen_lpair PROPERTY ENVIRONMENT "CEPGEN_PATH=${CEPGEN_PATH};PYTHONPATH=${CEPGEN_PATH}")
