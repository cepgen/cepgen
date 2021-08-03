if(NOT CMAKE_BUILD_UTILS)
  return()
endif()

message(STATUS "Utilities will be built")

#--- searching for LHAPDF
find_library(LHAPDF LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES lib)
find_path(LHAPDF_INCLUDE LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES include)

file(GLOB utils_standalone RELATIVE ${PROJECT_SOURCE_DIR}/src/utils *.cc)
if(NOT CMAKE_BUILD_CARD_CONV)
  file(GLOB cardconv RELATIVE ${PROJECT_SOURCE_DIR}/src/utils cepgenCardsConverter.cc)
  list(REMOVE_ITEM utils_standalone ${cardconv})
endif()
foreach(exec_src ${utils_standalone})
  get_filename_component(exec_bin ${exec_src} NAME_WE)
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES})
endforeach()

#-- build the utils requiring ROOT

if(ROOT_FOUND)
  foreach(exec_src ${utils_root})
    get_filename_component(exec_bin ${exec_src} NAME_WE)
    add_executable(${exec_bin} ${exec_src})
    set_property(TARGET ${exec_bin} PROPERTY CXX_STANDARD ${ROOT_CXX_STANDARD})
    target_include_directories(${exec_bin} PUBLIC ${ROOT_INCLUDE_DIRS})
    target_link_libraries(${exec_bin} PUBLIC ${CEPGEN_LIBRARIES} ${ROOT_LIBRARIES})
  endforeach()
endif()

