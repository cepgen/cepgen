cepgen_test_category(NAME "Physics"
                     SOURCES *.cc
                     FSOURCES *.f)

file(GLOB_RECURSE test_input_cards RELATIVE ${CMAKE_CURRENT_LIST_DIR} test_processes/* __init__.py)
foreach(_files ${test_input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${CMAKE_CURRENT_LIST_DIR}/test_processes.cfg test_processes.cfg COPYONLY)

add_test(NAME cross_section_vegas COMMAND processes_cross_sections -i MISER -q --cfg test_processes.cfg)
add_test(NAME cross_section_miser COMMAND processes_cross_sections -i Vegas -q --cfg test_processes.cfg)
set_tests_properties(cross_section_vegas cross_section_miser PROPERTIES WILL_FAIL TRUE)  # FIXME remove once agreement is reached
