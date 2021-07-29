file(GLOB cards_sources Handler.cpp LpairHandler.cpp CommandLineHandler.cpp)

include_directories(${PROJECT_SOURCE_DIR})

set(cards_libraries)
set(cards_definitions)
set(cards_includes)

#----- add Python parser if found

find_package(PythonInterp)
file(GLOB python_sources Python*.cpp)
if(CMAKE_VERSION VERSION_GREATER 3.12)
  find_package(Python COMPONENTS Interpreter Development)
  if(Python_FOUND)
    message(STATUS "Python v${Python_VERSION} found")
    list(APPEND cards_libraries ${Python_LIBRARIES})
    list(APPEND cards_definitions -DPYTHON)
    list(APPEND cards_sources ${python_sources})
    list(APPEND cards_includes ${Python_INCLUDE_DIRS})
  endif()
else()
  find_package(PythonLibs)
  if(PYTHONLIBS_FOUND)
    message(STATUS "Python v${PYTHONLIBS_VERSION_STRING} found")
    list(APPEND cards_libraries ${PYTHON_LIBRARIES})
    list(APPEND cards_definitions -DPYTHON)
    list(APPEND cards_sources ${python_sources})
    list(APPEND cards_includes ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
  endif()
endif()

#----- build the object

add_library(CepGenCards SHARED ${cards_sources})
target_compile_definitions(CepGenCards PUBLIC ${cards_definitions})
target_link_libraries(CepGenCards ${cards_libraries})
target_include_directories(CepGenCards PUBLIC ${cards_includes})

#----- installation rules

install(TARGETS CepGenCards DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lib)
