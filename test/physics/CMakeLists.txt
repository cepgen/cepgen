cepgen_test_category(NAME "Physics"
                     SOURCES *.cc
                     FSOURCES *.f)

file(GLOB_RECURSE test_input_cards RELATIVE ${CMAKE_CURRENT_LIST_DIR} test_processes/* __init__.py)
foreach(_files ${test_input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${CMAKE_CURRENT_LIST_DIR}/test_processes.cfg test_processes.cfg COPYONLY)

# one test per line of test_processes.cfg, so that they can be run in parallel (ctest -j) and individually (ctest -R)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CMAKE_CURRENT_LIST_DIR}/test_processes.cfg)
file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/test_processes.cfg _xsec_lines)
foreach(_line ${_xsec_lines})
  string(STRIP "${_line}" _line)
  if(_line STREQUAL "" OR _line MATCHES "^#")
    continue()
  endif()
  string(REGEX MATCH "^[^ \t]+" _name "${_line}")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/xsec_cfg/${_name}.cfg "${_line}\n")
  add_test(NAME cross_section_${_name}
           COMMAND processes_cross_sections -i Vegas -q --cfg xsec_cfg/${_name}.cfg)
  set_tests_properties(cross_section_${_name} PROPERTIES LABELS cross_sections)
endforeach()
