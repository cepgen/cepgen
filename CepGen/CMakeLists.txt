include_directories(${PROJECT_SOURCE_DIR})
set(PROC_PATH ${PROJECT_SOURCE_DIR}/CepGenProcesses)

file(GLOB core_sources Core/*.cpp Modules/*.cpp Utils/*.cpp *.cpp)
file(GLOB phys_sources Physics/*.cpp)
file(GLOB procdef_sources Processes/*.cpp Processes/Fortran/*.f)
file(GLOB sf_sources StructureFunctions/*.cpp)
file(GLOB proc_sources ${PROC_PATH}/*.cpp ${PROC_PATH}/*.f)

set(core_libraries ${GSL_LIB} ${GSL_CBLAS_LIB})
set(core_definitions)
set(core_includes)
set(proc_definitions)

#----- check the external dependencies for SFs/alphaS

if(GRV_SRC)
  message(STATUS "GRV PDFset found in ${GRV_SRC}!")
  list(APPEND sf_sources ${GRV_SRC})
  list(APPEND core_definitions -DGRVPDF)
endif()
if(LHAPDF)
  message(STATUS "LHAPDF found in ${LHAPDF}")
  list(APPEND core_libraries ${LHAPDF})
  list(APPEND core_includes ${LHAPDF_INCLUDE})
else()
  file(GLOB partonic_sf StructureFunctions/Partonic.cpp)
  list(REMOVE_ITEM sf_sources ${partonic_sf})
  file(GLOB lhapdf_alphas Physics/AlphaSLHAPDF.cpp)
  list(REMOVE_ITEM phys_sources ${lhapdf_alphas})
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
  list(APPEND core_libraries ${APFEL})
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
  list(APPEND core_definitions -DFUNC_MUPARSER)
  list(APPEND core_libraries ${MUPARSER})
elseif(EXPRTK)
  message(STATUS "Formulas parsing will use exprtk found in ${EXPRTK}")
  list(APPEND core_definitions -DFUNC_EXPRTK)
  list(APPEND core_includes ${EXPRTK}) # header-only
elseif(ROOT_FOUND)
  message(STATUS "Formulas parsing will use ROOT")
  list(APPEND core_definitions -DFUNC_ROOT)
  list(APPEND core_includes ${ROOT_INCLUDE_DIRS})
  list(APPEND core_libraries ${ROOT_LIBRARIES})
endif()

#----- build the objects

add_library(CepGenCore SHARED ${core_sources} ${phys_sources} ${procdef_sources} ${sf_sources})
target_link_libraries(CepGenCore ${core_libraries})
target_compile_definitions(CepGenCore PUBLIC ${core_definitions})
target_include_directories(CepGenCore PUBLIC ${core_includes} Processes/Fortran)

add_library(CepGenProcesses SHARED ${proc_sources})
target_compile_definitions(CepGenProcesses PUBLIC ${proc_definitions})
target_include_directories(CepGenProcesses PUBLIC Processes/Fortran)

#----- installation rules

install(TARGETS CepGenCore DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)
install(TARGETS CepGenProcesses DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)

