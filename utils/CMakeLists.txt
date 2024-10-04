if(NOT CMAKE_BUILD_UTILS)
  return()
endif()
message(STATUS "... utilities will be built")


#-- build all utils and link them to the core library
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")
file(GLOB sources RELATIVE ${CEPGEN_PATH}/utils *.cc)
foreach(exec_src ${sources})
  get_filename_component(exec_bin ${exec_src} NAME_WE)
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES})
endforeach()
