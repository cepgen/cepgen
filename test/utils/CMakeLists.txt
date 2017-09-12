#----- list all utilitaries

file(GLOB executables_noroot RELATIVE ${PROJECT_SOURCE_DIR}/test/utils *.cpp)
file(GLOB executables_root RELATIVE ${PROJECT_SOURCE_DIR}/test/utils *.cxx)

#----- build all utilitaries and link them to the core library

#----- start by the utils with no ROOT dependencies

foreach(exec_src ${executables_noroot})
  string(REPLACE ".cpp" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
  target_link_libraries(${exec_bin} ${CEPGEN_EXTERNAL_REQS} ${CEPGEN_LIBRARIES})
endforeach()

#----- specify the utilitaries requiring ROOT

find_package(ROOT QUIET)
if(ROOT_FOUND)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  endif()
  include_directories(${ROOT_INCLUDE_DIRS})
  link_directories(${ROOT_LIBRARY_DIR})

  foreach(exec_src ${executables_root})
    string(REPLACE ".cxx" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
    target_link_libraries(${exec_bin} ${CEPGEN_EXTERNAL_REQS} ${CEPGEN_LIBRARIES} ${ROOT_LIBRARIES})
  endforeach()
else()
  message(STATUS "ROOT not found! skipping these tests!")
endif()
