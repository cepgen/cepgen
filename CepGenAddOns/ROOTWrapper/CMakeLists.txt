#--- searching for ROOT
find_package(ROOT COMPONENTS Hist MathCore)
if(NOT ROOT_FOUND)
  return()
endif()
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  message(WARNING "ROOT bindings compilation is not yet supported...")
  return()
endif()

message(STATUS "... ROOT v${ROOT_VERSION} standard is c++${ROOT_CXX_STANDARD}")

set(sources IntegratorROOT.cpp
            ROOTDrawer.cpp
            ROOTTree*.cpp
            ROOTHistsHandler.cpp
            ROOTFunctional.cpp
            ROOTAnalyticalIntegrator.cpp)

string(REPLACE "-Qunused-arguments" "" ROOT_CXX_FLAGS ${ROOT_CXX_FLAGS})

#--- FOAM integrator
if(CMAKE_BUILD_FOAM)
  find_package(ROOT COMPONENTS Hist MathCore Foam)
  list(FIND ROOT_LIBRARIES "/usr/lib64/root/libFoam.so" _idx)
  if(${_idx} GREATER -1)
    list(APPEND sources IntegratorFoam.cpp)
  endif()
endif()

#----- build the object

cepgen_build(CepGenRoot SOURCES ${sources}
    EXT_HEADERS ${ROOT_INCLUDE_DIRS}
    EXT_LIBS ${ROOT_LIBRARIES}
    DEFINITIONS ${ROOT_DEFINITIONS}
    TESTS test/*.cc
    PROPERTY CXX_STANDARD ${ROOT_CXX_STANDARD}
    CFLAGS "${ROOT_CXX_FLAGS} -Wno-bitwise-instead-of-logical -Wno-sizeof-pointer-div -Wno-unused-command-line-argument"
    INSTALL_COMPONENT root)
