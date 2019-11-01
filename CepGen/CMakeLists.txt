file(GLOB core_sources Core/*.cpp Modules/*.cpp *.cpp)
file(GLOB phys_sources Physics/*.cpp)
file(GLOB sf_sources StructureFunctions/*.cpp)

include_directories(${PROJECT_SOURCE_DIR})

message(STATUS "GSL found in ${GSL_LIB}")
message(STATUS "GSL CBLAS found in ${GSL_CBLAS_LIB}")

set(core_libraries ${GSL_LIB} ${GSL_CBLAS_LIB})
set(core_definitions)
set(core_includes)

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

#--- linking with LHAPDF

if(LHAPDF)
  message(STATUS "LHAPDF found in ${LHAPDF}")
  list(APPEND core_libraries ${LHAPDF})
  list(APPEND core_includes ${LHAPDF_INCLUDE})
else()
  file(GLOB partonic_sf StructureFunctions/Partonic.cpp)
  list(REMOVE_ITEM sf_sources ${partonic_sf})
endif()

#----- build the objects

add_library(CepGenCore SHARED ${core_sources} ${phys_sources} ${sf_sources})
target_link_libraries(CepGenCore ${core_libraries})
target_compile_definitions(CepGenCore PUBLIC ${core_definitions})
target_include_directories(CepGenCore PUBLIC ${core_includes})

#----- installation rules

install(TARGETS CepGenCore DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT lib)

