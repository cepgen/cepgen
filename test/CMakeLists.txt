#-- list all test executables

file(GLOB tests_standalone RELATIVE ${PROJECT_SOURCE_DIR}/test *.cc)
file(GLOB tests_fortran RELATIVE ${PROJECT_SOURCE_DIR}/test *.f)

#-- include the utilitaries

set(tests)
set(tests_libraries ${CEPGEN_LIBRARIES})

#-- skip functionals tests if no parsing library is found

if(NOT CEPGEN_HAS_FUNCTIONAL)
  message(STATUS "No functionals parser found on this system. Disabling the corresponding tests.")
  file(GLOB tests_functionals
      RELATIVE ${PROJECT_SOURCE_DIR}/test
      test_integrator.cc test_function_parser*.*)
#list(REMOVE_ITEM tests_standalone ${tests_functionals})
endif()

#-- build all tests and link them to the core library

foreach(exec_src ${tests_standalone})
  get_filename_component(exec_bin ${exec_src} NAME_WE)
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${tests_libraries})
  if(CEPGEN_ROOT_FUNCTIONALS)
    set_property(TARGET ${exec_bin} PROPERTY CXX_STANDARD ${ROOT_CXX_STANDARD})
  endif()
  list(APPEND tests ${exec_bin})
endforeach()

foreach(exec_src ${tests_fortran})
  get_filename_component(exec_bin ${exec_src} NAME_WE)
  add_executable(${exec_bin} ${exec_src})
  target_link_libraries(${exec_bin} ${tests_libraries} stdc++)
#list(APPEND tests ${exec_bin}) #for the moment, excluding F77 tests
endforeach()

file(GLOB_RECURSE test_input_cards RELATIVE ${PROJECT_SOURCE_DIR}/test test_processes/* __init__.py)
foreach(_files ${test_input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${PROJECT_SOURCE_DIR}/test/test_processes.cfg test_processes.cfg COPYONLY)

foreach(test ${tests})
  message(STATUS "Added test ${test}")
  add_test(NAME ${test} COMMAND ${test} -h)
endforeach()
