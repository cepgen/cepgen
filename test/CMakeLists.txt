#----- list all test executables

file(GLOB tests_noroot RELATIVE ${PROJECT_SOURCE_DIR}/test *.cc)
file(GLOB tests_root RELATIVE ${PROJECT_SOURCE_DIR}/test *.cxx)
file(GLOB tests_fortran RELATIVE ${PROJECT_SOURCE_DIR}/test *.f)

#----- include the utilitaries

set(tests)
set(tests_libraries ${CEPGEN_LIBRARIES})

file(GLOB tests_utils_sources ArgumentsParser.cpp)
add_library(CepGenTestUtils SHARED ${tests_utils_sources})
list(APPEND tests_libraries CepGenTestUtils)

#----- build all tests and link them to the core library

foreach(exec_src ${tests_noroot})
  string(REPLACE ".cc" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${tests_libraries})
  list(APPEND tests ${exec_bin})
endforeach()

foreach(exec_src ${tests_fortran})
  string(REPLACE ".f" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${tests_libraries})
  #list(APPEND tests ${exec_bin}) # for the moment, excluding F77 tests
endforeach()

#----- specify the tests requiring ROOT

if(ROOT_FOUND)
  include_directories(${ROOT_INCLUDE_DIRS})
  link_directories(${ROOT_LIBRARY_DIR})

  foreach(exec_src ${tests_root})
    string(REPLACE ".cxx" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_property(TARGET ${exec_bin} PROPERTY CXX_STANDARD 14)
    target_link_libraries(${exec_bin} ${tests_libraries} ${ROOT_LIBRARIES})
    list(APPEND tests ${exec_bin})
  endforeach()
else()
  message(STATUS "ROOT not found! skipping these tests!")
endif()

file(GLOB_RECURSE test_input_cards RELATIVE ${PROJECT_SOURCE_DIR}/test test_processes/* __init__.py)
foreach(_files ${test_input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${PROJECT_SOURCE_DIR}/test/test_processes.cfg test_processes.cfg COPYONLY)

foreach(test ${tests})
  if(${CMAKE_BUILD_TESTS})
    set_target_properties(${test} PROPERTIES EXCLUDE_FROM_ALL false)
  else()
    set_target_properties(${test} PROPERTIES EXCLUDE_FROM_ALL true)
  endif()
  add_test(NAME ${test} COMMAND ${test} -h)
  message(STATUS "Added test ${test}")
endforeach()
