find_package(Boost)
if(NOT Boost_FOUND)
  return()
endif()

set(sources BoostTreeHandler.cpp
            BoostTreeUtils.cpp
            BoostTrapAnalyticalIntegrator.cpp
            BoostGaussKronrodAnalyticalIntegrator.cpp
            BoostGaussLegendreAnalyticalIntegrator.cpp)
file(GLOB boost_int_include ${Boost_INCLUDE_DIRS}/boost/math/quadrature/naive_monte_carlo.hpp)
if(boost_int_include)
  list(APPEND sources IntegratorNaive.cpp)
endif()

#----- build the object

cepgen_build(CepGenBoost SOURCES ${sources}
    EXT_LIBS ${Boost_LIBRARIES}
    EXT_HEADERS ${Boost_INCLUDE_DIRS}
    CFLAGS "-pthread -Wno-pedantic"
    INSTALL_COMPONENT boost)
cpack_add_component(boost
    DISPLAY_NAME "CepGen Boost wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Boost library"
    DEPENDS lib)

#----- Python bindings

find_package(Boost COMPONENTS system python)
find_package(Python COMPONENTS Interpreter Development)
if(Python_FOUND)
  add_library(pycepgen SHARED PythonWrapper.cpp PythonUtils.cpp)
  target_link_libraries(pycepgen CepGen ${Boost_LIBRARIES})
  target_include_directories(pycepgen PUBLIC ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
  file(TOUCH_NOCREATE ${CMAKE_CURRENT_BINARY_DIR}/src/__init__.py)
  set_target_properties(pycepgen PROPERTIES PREFIX "")
  target_include_directories(pycepgen PUBLIC ${Python_INCLUDE_DIRS})
  install(TARGETS pycepgen DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT boost)
endif()
