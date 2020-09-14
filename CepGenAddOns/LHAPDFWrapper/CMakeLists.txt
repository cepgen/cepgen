#--- searching for LHAPDF
find_library(LHAPDF LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES lib)
find_path(LHAPDF_INCLUDE LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES include)

if(LHAPDF)
  message(STATUS "LHAPDF found in ${LHAPDF}")
else()
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources PartonicStructureFunctionsLHAPDF.cpp AlphaSLHAPDF.cpp)

#----- build the objects

add_library(CepGenLHAPDF SHARED ${sources})
target_link_libraries(CepGenLHAPDF PRIVATE ${LHAPDF})
target_include_directories(CepGenLHAPDF PRIVATE ${LHAPDF_INCLUDE})
list(APPEND CEPGEN_ADDONS_EXT CepGenLHAPDF)
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenLHAPDF DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lhapdf)
