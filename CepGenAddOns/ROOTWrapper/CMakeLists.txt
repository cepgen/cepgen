#--- searching for ROOT
find_package(ROOT COMPONENTS Hist MathCore)
if(NOT ROOT_FOUND)
  return()
endif()
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  message(WARNING "ROOT bindings compilation is not yet supported...")
  return()
endif()

file(GLOB root_sources ROOT*.cpp
                       AnalyticalIntegrator.cpp
                       Derivator.cpp
                       Drawer.cpp
                       EventExporter.cpp
                       EventImporter.cpp
                       Functional.cpp
                       Integrator.cpp
                       RandomGenerator.cpp)

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
cpack_add_component(root
    DISPLAY_NAME "CepGen ROOT wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the ROOT library"
    DEPENDS lib)

message(STATUS "... ROOT v${ROOT_VERSION} standard is c++${ROOT_CXX_STANDARD}")
