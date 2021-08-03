#--- searching for ROOT
find_package(ROOT COMPONENTS Hist MathCore Foam)
if(NOT ROOT_FOUND)
  return()
endif()

set(sources IntegratorROOT.cpp ROOTTree*.cpp ROOTHistsHandler.cpp FunctionalROOT.cpp)

#--- FOAM integrator
if(CMAKE_BUILD_FOAM)
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
    PROPERTY CXX_STANDARD ${ROOT_CXX_STANDARD}
    CFLAGS ${ROOT_CXX_FLAGS}
    INSTALL_COMPONENT root)
message(STATUS "... ROOT v${ROOT_VERSION} standard is c++${ROOT_CXX_STANDARD}")

message(STATUS "Formulas parsing will use ROOT found in ${ROOT_INCLUDE_DIRS}")
set(CEPGEN_HAS_FUNCTIONAL ON PARENT_SCOPE)

