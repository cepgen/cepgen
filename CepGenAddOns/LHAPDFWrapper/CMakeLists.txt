if(LHAPDF)
  message(STATUS "LHAPDF found in ${LHAPDF}")
else()
  message(FATAL_ERROR "Trying to build LHAPDF wrapper library while no LHAPDF installation was found on system!")
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
