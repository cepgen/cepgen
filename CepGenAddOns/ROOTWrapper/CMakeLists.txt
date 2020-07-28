if(NOT ROOT_FOUND)
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB int_sources IntegratorROOT.cpp)
file(GLOB out_sources ROOTTree*.cpp ROOTHistsHandler.cpp)
file(GLOB utils_sources FunctionalROOT.cpp)

list(FIND ROOT_LIBRARIES "/usr/lib64/root/libFoam.so" _idx)
if(${_idx} GREATER -1)
  list(APPEND int_sources IntegratorFoam.cpp)
endif()

#----- build the objects

add_library(CepGenRoot SHARED ${int_sources} ${out_sources} ${utils_sources})
target_compile_definitions(CepGenRoot PRIVATE ${ROOT_DEFINITIONS})
target_link_libraries(CepGenRoot PRIVATE ${ROOT_LIBRARIES})
target_include_directories(CepGenRoot PRIVATE ${ROOT_INCLUDE_DIRS})
set_property(TARGET CepGenRoot PROPERTY CXX_STANDARD 14)
list(APPEND CEPGEN_ADDONS_EXT CepGenRoot)
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

message(STATUS "Formulas parsing will use ROOT found in ${ROOT_INCLUDE_DIRS}")
set(CEPGEN_HAS_FUNCTIONALS ON PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenRoot DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT root)
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT root FILES_MATCHING PATTERN "*.h")
