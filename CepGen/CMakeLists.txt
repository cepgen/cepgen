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
        CollinearFluxes/*.cpp
        Event/*.cpp
        FormFactors/*.cpp
        KTFluxes/*.cpp
        Physics/*.cpp
        StructureFunctions/*.cpp)
file(GLOB proc_sources
        Process/*.cpp
        Process/Fortran/*.f)
file(GLOB core_includes
        Cards
        CollinearFluxes
        Core
        EventFilter
        FormFactors
        Integration
        KTFluxes
        Modules
        OutputModules
        Physics
        Process
        StructureFunctions
        Utils)

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

include(FindVersion)
add_version_definition(${core_sources})

cepgen_build(CepGen
        SOURCES ${core_sources} ${phys_sources} ${proc_sources}
        EXT_LIBS ${CEPGEN_CORE_EXT}
        EXT_HEADERS ${core_includes}
        INSTALL_COMPONENT lib)
