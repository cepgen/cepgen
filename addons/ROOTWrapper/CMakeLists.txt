#--- searching for ROOT
find_package(ROOT COMPONENTS Hist MathCore OPTIONAL_COMPONENTS Foam)
if(NOT ROOT_FOUND)
  return()
endif()
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  message(WARNING "ROOT bindings compilation is not yet supported...")
  return()
endif()

file(GLOB root_sources src/ROOT*.cpp
                       src/AnalyticalIntegrator.cpp
                       src/Derivator.cpp
                       src/Drawer.cpp
                       src/EventExporter.cpp
                       src/EventImporter.cpp
                       src/Functional.cpp
                       src/Integrator.cpp
                       src/RandomGenerator.cpp)

#--- FOAM integrator
if("${ROOT_LIBRARIES}" MATCHES "^?(Foam)?$")
  file(GLOB foam_sources src/Foam*.cpp)
endif()

#----- build the object
cepgen_build(CepGenRoot
    SOURCES ${root_sources} ${foam_sources}
    LIBRARIES ${ROOT_LIBRARIES}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS}
    TESTS test/*.cc
    PROPERTY CXX_STANDARD ${ROOT_CXX_STANDARD}
    CFLAGS "${ROOT_CXX_FLAGS}"
    COMPONENT root)
cpack_add_component(root
    DISPLAY_NAME "CepGen ROOT wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the ROOT library"
    DEPENDS lib)
install(DIRECTORY CepGenRoot
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT root)

message(STATUS "... ROOT v${ROOT_VERSION} standard is c++${ROOT_CXX_STANDARD}")
