#----- list all test executables

file(GLOB executables_noroot RELATIVE ${PROJECT_SOURCE_DIR}/test *.cpp)
file(GLOB executables_root RELATIVE ${PROJECT_SOURCE_DIR}/test *.cxx)
file(GLOB executables_fortran RELATIVE ${PROJECT_SOURCE_DIR}/test *.f)

#----- include the utilitaries

set(tests "")

#----- build all tests and link them to the core library

foreach(exec_src ${executables_noroot})
  string(REPLACE ".cpp" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
  target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES})
  list(APPEND tests ${exec_bin})
endforeach()

foreach(exec_src ${executables_fortran})
  string(REPLACE ".f" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
  target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES})
  list(APPEND tests ${exec_bin})
endforeach()

#----- specify the tests requiring ROOT

find_package(ROOT QUIET)
if(ROOT_FOUND)
  set(CMAKE_CXX_STANDARD 14)
  include_directories(${ROOT_INCLUDE_DIRS})
  include_directories(${ROOT_INCLUDE_DIR})
  link_directories(${ROOT_LIBRARY_DIR})

  foreach(exec_src ${executables_root})
    string(REPLACE ".cxx" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
    target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES} ${ROOT_LIBRARIES})
    list(APPEND tests ${exec_bin})
  endforeach()
else()
  message(STATUS "ROOT not found! skipping these tests!")
endif()

message(STATUS "list of available tests: ${tests}")
