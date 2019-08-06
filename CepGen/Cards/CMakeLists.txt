file(GLOB cards_sources Handler.cpp LpairHandler.cpp)

include_directories(${PROJECT_SOURCE_DIR})

set(cards_libraries "")

#----- add Python parser if found

if(PYTHONLIBS_FOUND)
  list(APPEND cards_libraries ${PYTHON_LIBRARIES})
  file(GLOB python_sources Python*.cpp)
  add_definitions(-DPYTHON)
  list(APPEND cards_sources ${python_sources})
  message(STATUS "Python v${PYTHONLIBS_VERSION_STRING} found")
  include_directories(${PYTHON_INCLUDE_DIRS})
endif()

#----- build the object

add_library(CepGenCards SHARED ${cards_sources})
target_link_libraries(CepGenCards ${cards_libraries})

#----- installation rules

install(TARGETS CepGenCards DESTINATION ${LIB_INSTALL_DIR} COMPONENT lib)

