if(NOT CMAKE_BUILD_UTILS)
  return()
endif()

message(STATUS "Utilities will be built")

#-- build all utils and link them to the core library

file(GLOB utils_standalone RELATIVE ${PROJECT_SOURCE_DIR}/src/utils *.cc)
foreach(exec_src ${utils_standalone})
  get_filename_component(exec_bin ${exec_src} NAME_WE)
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES})
endforeach()

#-- specify the utils requiring ROOT

if(ROOT_FOUND)
  file(GLOB utils_root RELATIVE ${PROJECT_SOURCE_DIR}/src/utils *.cxx)
  foreach(exec_src ${utils_root})
    get_filename_component(exec_bin ${exec_src} NAME_WE)
    add_executable(${exec_bin} ${exec_src})
    set_property(TARGET ${exec_bin} PROPERTY CXX_STANDARD ${ROOT_CXX_STANDARD})
    target_include_directories(${exec_bin} PUBLIC ${ROOT_INCLUDE_DIRS})
    target_link_libraries(${exec_bin} PUBLIC ${CEPGEN_LIBRARIES} ${ROOT_LIBRARIES})
  endforeach()
endif()

