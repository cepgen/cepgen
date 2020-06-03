include_directories(${PROJECT_SOURCE_DIR})

file(GLOB core_sources Core/*.cpp Modules/*.cpp Utils/*.cpp *.cpp)
file(GLOB phys_sources Event/*.cpp Physics/*.cpp)
file(GLOB int_sources Integration/*.cpp)
file(GLOB proc_sources Processes/*.cpp Processes/Fortran/*.f)
file(GLOB sf_sources StructureFunctions/*.cpp)

set(core_definitions)
set(core_includes)

#----- check the external dependencies for SFs/alphaS

if(GRV_SRC)
  message(STATUS "GRV PDFset found in ${GRV_SRC}!")
  list(APPEND sf_sources ${GRV_SRC})
else()
  message(STATUS "GRV PDFset not found. Will not build Szczurek-Uleshchenko SFs")
  file(GLOB su_sf StructureFunctions/SzczurekUleshchenko.cpp)
  list(REMOVE_ITEM sf_sources ${su_sf})
endif()
if(ALPHAS_SRC)
  message(STATUS "alphaS evolution algorithm found in ${ALPHAS_SRC}")
  list(APPEND phys_sources ${ALPHAS_SRC})
else()
  file(GLOB pegasus_alphas Physics/AlphaSPEGASUS.cpp)
  file(REMOVE_ITEM phys_sources ${pegasus_alphas})
endif()
if(APFEL)
  message(STATUS "APFEL found in ${APFEL}")
  list(APPEND core_includes ${APFEL_INCLUDE})
  list(APPEND CEPGEN_CORE_EXT ${APFEL})
else()
  file(GLOB apfel_alphas Physics/AlphaSAPFEL.cpp)
  list(REMOVE_ITEM phys_sources ${apfel_alphas})
endif()

#----- check the external dependencies for physics utilities

if(KMR_INTERP)
  message(STATUS "KMR interpolation algorithm found in ${KMR_INTERP}")
  list(APPEND phys_sources ${KMR_INTERP})
endif()

if(MUPARSER)
  message(STATUS "Formulas parsing will use muParser found in ${MUPARSER}")
  list(APPEND CEPGEN_CORE_EXT ${MUPARSER})
  set(CEPGEN_HAS_FUNCTIONALS ON)
else()
  file(GLOB func_muparser Utils/FunctionalMuParser.cpp)
  list(REMOVE_ITEM core_sources ${func_muparser})
endif()

if(EXPRTK)
  message(STATUS "Formulas parsing will use exprtk found in ${EXPRTK}")
  list(APPEND core_includes ${EXPRTK}) # header-only
  set(CEPGEN_HAS_FUNCTIONALS ON)
else()
  file(GLOB func_exprtk Utils/FunctionalExprTk.cpp)
  list(REMOVE_ITEM core_sources ${func_exprtk})
endif()

set(CEPGEN_HAS_FUNCTIONALS ${CEPGEN_HAS_FUNCTIONALS} PARENT_SCOPE)

#----- build the objects

add_library(CepGenCore SHARED ${core_sources} ${int_sources} ${phys_sources} ${proc_sources} ${sf_sources})
target_compile_definitions(CepGenCore PUBLIC ${core_definitions})
target_link_libraries(CepGenCore PUBLIC ${CEPGEN_CORE_EXT})
target_include_directories(CepGenCore PUBLIC ${core_includes})

#----- installation rules

install(TARGETS CepGenCore DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lib)
