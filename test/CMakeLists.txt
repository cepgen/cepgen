#----- list all test executables

file(GLOB executables_noroot RELATIVE ${PROJECT_SOURCE_DIR}/test *.cpp)
file(GLOB executables_root RELATIVE ${PROJECT_SOURCE_DIR}/test *.cxx)
file(GLOB executables_fortran RELATIVE ${PROJECT_SOURCE_DIR}/test *.f)

#----- include the utilitaries

set(tests)
set(tests_libraries ${CEPGEN_LIBRARIES})

file(GLOB tests_utils_sources ArgumentsParser.cpp)
add_library(CepGenTestUtils OBJECT ${tests_utils_sources})
list(APPEND tests_libraries CepGenTestUtils)

#----- build all tests and link them to the core library

foreach(exec_src ${executables_noroot})
  string(REPLACE ".cpp" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
  target_link_libraries(${exec_bin} ${tests_libraries})
  list(APPEND tests ${exec_bin})
endforeach()

foreach(exec_src ${executables_fortran})
  string(REPLACE ".f" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
  target_link_libraries(${exec_bin} ${tests_libraries})
  list(APPEND tests ${exec_bin})
endforeach()

#----- specify the tests requiring ROOT

if(ROOT_FOUND)
  include_directories(${ROOT_INCLUDE_DIRS})
  link_directories(${ROOT_LIBRARY_DIR})

  foreach(exec_src ${executables_root})
    string(REPLACE ".cxx" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
    set_property(TARGET ${exec_bin} PROPERTY CXX_STANDARD 14)
    target_link_libraries(${exec_bin} ${tests_libraries} ${ROOT_LIBRARIES})
    list(APPEND tests ${exec_bin})
  endforeach()
else()
  message(STATUS "ROOT not found! skipping these tests!")
endif()

message(STATUS "The following tests/executables are available:\n${tests}")

file(GLOB_RECURSE test_input_cards RELATIVE ${PROJECT_SOURCE_DIR}/test test_processes/* __init__.py)
foreach(_files ${test_input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${PROJECT_SOURCE_DIR}/test/test_processes.cfg test_processes.cfg COPYONLY)
