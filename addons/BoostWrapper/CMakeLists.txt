find_package(Boost)
if(NOT Boost_FOUND)
  return()
endif()

set(sources src/BoostTreeHandler.cpp
            src/BoostTreeUtils.cpp
            src/FiniteDifferenceDerivator.cpp
            src/TrapIntegrator.cpp
            src/GaussKronrodIntegrator.cpp
            src/GaussLegendreIntegrator.cpp)
file(GLOB boost_int_include ${Boost_INCLUDE_DIRS}/boost/math/quadrature/naive_monte_carlo.hpp)
if(boost_int_include)
  list(APPEND sources src/NaiveIntegrator.cpp)
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
