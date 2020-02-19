#----- list all test executables

file(GLOB tests_standalone RELATIVE ${PROJECT_SOURCE_DIR}/test *.cc)
file(GLOB tests_root RELATIVE ${PROJECT_SOURCE_DIR}/test *.cxx)
file(GLOB tests_fortran RELATIVE ${PROJECT_SOURCE_DIR}/test *.f)
file(GLOB tests_functionals RELATIVE ${PROJECT_SOURCE_DIR}/test test_function_parser*.*)

#----- include the utilitaries

set(tests)
set(tests_libraries ${CEPGEN_LIBRARIES})

#----- skip functionals tests if no parsing library is found

if(NOT CEPGEN_HAS_FUNCTIONALS)
  list(REMOVE_ITEM tests_standalone ${tests_functionals})
  list(REMOVE_ITEM tests_root ${tests_functionals})
endif()

#----- build all tests and link them to the core library

foreach(exec_src ${tests_standalone})
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
  foreach(exec_src ${tests_root})
    string(REPLACE ".cxx" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_property(TARGET ${exec_bin} PROPERTY CXX_STANDARD 14)
    target_include_directories(${exec_bin} PUBLIC ${ROOT_INCLUDE_DIRS})
    target_link_libraries(${exec_bin} PUBLIC ${tests_libraries} ${ROOT_LIBRARIES})
    list(APPEND tests ${exec_bin})
  endforeach()
endif()

file(GLOB_RECURSE test_input_cards RELATIVE ${PROJECT_SOURCE_DIR}/test test_processes/* __init__.py)
foreach(_files ${test_input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${PROJECT_SOURCE_DIR}/test/test_processes.cfg test_processes.cfg COPYONLY)

foreach(test ${tests})
  message(STATUS "Added test ${test}")
  add_test(NAME ${test} COMMAND ${test} -h)
endforeach()
