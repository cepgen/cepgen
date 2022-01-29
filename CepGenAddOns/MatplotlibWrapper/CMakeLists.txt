#--- searching for MATPLOTLIB
find_package(matplotlib_cpp)

if(NOT ${matplotlib_cpp_LIBRARIES})
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources DrawerMatplotlib.cpp MatplotlibHistsHandler.cpp)

#----- build the object

cepgen_build(CepGenMatplotlib SOURCES ${sources}
    EXT_LIBS ${matplotlib_cpp_LIBRARIES} matplotlib_cpp::matplotlib_cpp
    EXT_HEADERS ${matplotlib_cpp_INCLUDE_DIRS}
    INSTALL_COMPONENT matplotlib)
