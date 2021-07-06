#-- build all utils and link them to the core library

file(GLOB utils_standalone *.cc)
foreach(exec_src ${utils_standalone})
  string(REPLACE ".cc" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES})
endforeach()

#-- specify the utils requiring ROOT

find_package(ROOT QUIET)
if(ROOT_FOUND)
  file(GLOB utils_root RELATIVE ${. *.cxx)
  foreach(exec_src ${utils_root})
    string(REPLACE ".cxx" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_property(TARGET ${exec_bin} PROPERTY CXX_STANDARD 17)
    target_include_directories(${exec_bin} PUBLIC ${ROOT_INCLUDE_DIRS})
    target_link_libraries(${exec_bin} PUBLIC ${CEPGEN_LIBRARIES} ${ROOT_LIBRARIES})
  endforeach()
endif()

