#-- list all test executables

file(GLOB tests_standalone RELATIVE ${CMAKE_CURRENT_LIST_DIR} *.cc)
file(GLOB tests_fortran RELATIVE ${CMAKE_CURRENT_LIST_DIR} *.f)

#-- include the utilitaries

set(tests)
set(tests_libraries ${CEPGEN_LIBRARIES})

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

message(STATUS ">>> Utilities tests <<<")
foreach(test ${tests})
  message(STATUS "... added test ${test}")
  add_test(NAME ${test} COMMAND ${test} -h)
endforeach()
