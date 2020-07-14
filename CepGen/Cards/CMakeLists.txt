file(GLOB cards_sources Handler.cpp LpairHandler.cpp CommandLineHandler.cpp)

include_directories(${PROJECT_SOURCE_DIR})

set(cards_libraries)
set(cards_definitions)
set(cards_includes)

#----- add Python parser if found

if(PYTHONLIBS_FOUND)
  message(STATUS "Python v${PYTHONLIBS_VERSION_STRING} found")
  file(GLOB python_sources Python*.cpp)
  list(APPEND cards_libraries ${PYTHON_LIBRARIES})
  list(APPEND cards_definitions -DPYTHON)
  list(APPEND cards_sources ${python_sources})
  list(APPEND cards_includes ${PYTHON_INCLUDE_DIRS})
endif()

#----- build the object

add_library(CepGenCards SHARED ${cards_sources})
target_compile_definitions(CepGenCards PUBLIC ${cards_definitions})
target_link_libraries(CepGenCards ${cards_libraries})
target_include_directories(CepGenCards PUBLIC ${cards_includes})

#----- installation rules

install(TARGETS CepGenCards DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lib)
