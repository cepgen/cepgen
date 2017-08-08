#----- list all test executables

file(GLOB executables_noroot RELATIVE ${PROJECT_SOURCE_DIR}/test *.cpp)
file(GLOB executables_root RELATIVE ${PROJECT_SOURCE_DIR}/test *.cxx)

#----- build all tests and link them to the core library

foreach(exec_src ${executables_noroot})
    string(REPLACE ".cpp" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
    target_link_libraries(${exec_bin} ${GSL_LIB} ${CEPGEN_LIBRARIES})
    add_test(${exec_bin} ${exec_bin})
    set_tests_properties(${exec_bin} PROPERTIES PASS_REGULAR_EXPRESSION "Passed")
endforeach()

#----- specify the tests requiring ROOT

find_package(ROOT QUIET)
if(${ROOT_FOUND})
    message(STATUS "ROOT found in ${ROOT_INCLUDE_DIRS}")
    include_directories(${ROOT_INCLUDE_DIRS})
    link_directories(${ROOT_LIBRARY_DIR})

    foreach(exec_src ${executables_root})
        string(REPLACE ".cxx" "" exec_bin ${exec_src})
        add_executable(${exec_bin} ${exec_src})
        set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
        target_link_libraries(${exec_bin} ${GSL_LIB} ${CEPGEN_LIBRARIES} ${ROOT_LIBRARIES})
    endforeach()
else()
    message(STATUS "ROOT not found! skipping these tests!")
endif()
