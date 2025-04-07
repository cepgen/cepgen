find_package(Boost)
if(NOT Boost_FOUND)
  return()
endif()

set(sources src/BoostTreeHandler.cpp
            src/BoostTreeUtils.cpp
            src/BoostTrapAnalyticalIntegrator.cpp
            src/BoostGaussKronrodAnalyticalIntegrator.cpp
            src/BoostGaussLegendreAnalyticalIntegrator.cpp
            src/FiniteDifferenceDerivator.cpp)
file(GLOB boost_int_include ${Boost_INCLUDE_DIRS}/boost/math/quadrature/naive_monte_carlo.hpp)
if(boost_int_include)
  list(APPEND sources src/IntegratorNaive.cpp)
endif()

#----- build the object

cepgen_build(CepGenBoost
    SOURCES ${sources}
    LIBRARIES ${Boost_LIBRARIES}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS}
    CFLAGS "-pthread -Wno-pedantic"
    COMPONENT boost)
cpack_add_component(boost
    DISPLAY_NAME "CepGen Boost wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Boost library"
    DEPENDS lib)
install(DIRECTORY CepGenBoost
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT boost)

#----- Python bindings

find_package(Boost COMPONENTS system python)
find_package(Python COMPONENTS Interpreter Development)
if(Python_FOUND)
  add_library(pycepgen SHARED src/PythonWrapper.cpp src/PythonUtils.cpp)
  target_link_libraries(pycepgen CepGen ${Boost_LIBRARIES})
  target_include_directories(pycepgen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
  file(TOUCH_NOCREATE ${CMAKE_CURRENT_BINARY_DIR}/src/__init__.py)
  set_target_properties(pycepgen PROPERTIES PREFIX "")
  target_include_directories(pycepgen PUBLIC ${Python_INCLUDE_DIRS})
  install(TARGETS pycepgen DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT boost)
endif()
