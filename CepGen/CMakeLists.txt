file(GLOB core_sources
        Core/*.cpp
        Modules/*.cpp
        Utils/*.cpp
        Integration/*.cpp
        *.cpp)
file(GLOB phys_sources
        Event/*.cpp
        Physics/*.cpp
        StructureFunctions/*.cpp
        FormFactors/*.cpp)
file(GLOB proc_sources
        Processes/*.cpp
        Processes/Fortran/*.f)
file(GLOB core_includes
        Core
        Modules
        Utils
        Physics
        Integration
        Processes
        StructureFunctions
        FormFactors)

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

cepgen_build(CepGenCore
        SOURCES ${core_sources} ${phys_sources} ${proc_sources}
        EXT_LIBS ${CEPGEN_CORE_EXT} stdc++fs
        EXT_HEADERS ${core_includes}
        INSTALL_COMPONENT lib)
