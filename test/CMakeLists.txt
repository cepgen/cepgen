#-- list all test executables

file(GLOB tests_standalone RELATIVE ${PROJECT_SOURCE_DIR}/test *.cc)
file(GLOB tests_fortran RELATIVE ${PROJECT_SOURCE_DIR}/test *.f)

#-- include the utilitaries

#--- searching for ROOT
find_package(ROOT QUIET)
#--- searching for LHAPDF
find_library(LHAPDF LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES lib)
find_path(LHAPDF_INCLUDE LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES include)

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

if(tests_lhapdf AND LHAPDF)
  list(REMOVE_ITEM tests_standalone ${tests_lhapdf})
  list(REMOVE_ITEM tests_root ${tests_lhapdf})
  foreach(exec_src ${tests_lhapdf})
    string(REPLACE ".cxx" "" exec_bin ${exec_src})
    string(REPLACE ".cc" "" exec_src ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_property(TARGET ${exec_bin} PROPERTY CXX_STANDARD 17)
    target_include_directories(${exec_bin} PUBLIC ${LHAPDF_INCLUDE} ${ROOT_INCLUDE_DIRS})
    target_link_libraries(${exec_bin} PUBLIC ${tests_libraries} ${LHAPDF} ${ROOT_LIBRARIES})
    list(APPEND tests ${exec_bin})
  endforeach()
endif()

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
  #list(APPEND tests ${exec_bin}) # for the moment, excluding F77 tests
endforeach()

message(STATUS ">>> Miscellaneous tests <<<")
foreach(test ${tests})
  message(STATUS "... added test ${test}")
  add_test(NAME ${test} COMMAND ${test} -h)
endforeach()
add_test(NAME integrator COMMAND test_integrator -q -a)

add_subdirectory(utils)
add_subdirectory(physics)
