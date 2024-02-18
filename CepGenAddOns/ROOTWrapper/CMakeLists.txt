#--- searching for ROOT
find_package(ROOT COMPONENTS Hist MathCore)
if(NOT ROOT_FOUND)
  return()
endif()
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  message(WARNING "ROOT bindings compilation is not yet supported...")
  return()
endif()

include(${ROOT_USE_FILE})
set(CMAKE_ROOTTEST_NOROOTMAP OFF)
file(GLOB EVENT_HEADERS ${PROJECT_SOURCE_DIR}/CepGen/Event/*.h)
set(ROOT_DICT_NAME "CepGenEventDict")
reflex_generate_dictionary(${ROOT_DICT_NAME} classes.h ${EVENT_HEADERS} SELECTION classes_def.xml)

file(GLOB root_sources ROOT*.cpp ${ROOT_DICT_NAME})
#add_dependencies(CepGenEvent CepGenEventDict-dictgen)

string(REPLACE "-Qunused-arguments" "" ROOT_CXX_FLAGS ${ROOT_CXX_FLAGS})

#--- FOAM integrator
if(CMAKE_BUILD_FOAM)
  find_package(ROOT COMPONENTS Hist MathCore Foam)
  list(FIND ROOT_LIBRARIES "/usr/lib64/root/libFoam.so" _idx)
  if(${_idx} GREATER -1)
    file(GLOB foam_sources Foam*.cpp)
  endif()
endif()

#----- build the object

cepgen_build(CepGenRoot SOURCES ${root_sources} ${foam_sources}
    EXT_HEADERS ${ROOT_INCLUDE_DIRS}
    EXT_LIBS ${ROOT_LIBRARIES}
    TESTS test/*.cc
    PROPERTY CXX_STANDARD ${ROOT_CXX_STANDARD}
    CFLAGS "${ROOT_CXX_FLAGS}"
    INSTALL_COMPONENT root)

message(STATUS "... ROOT v${ROOT_VERSION} standard is c++${ROOT_CXX_STANDARD}")
#add_dependencies(CepGenRoot ${ROOT_DICT_NAME}-dictgen)
