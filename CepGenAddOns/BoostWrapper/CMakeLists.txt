find_package(Boost)

if(NOT Boost_FOUND)
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

set(sources BoostTreeHandler.cpp)
file(GLOB boost_int_include ${Boost_INCLUDE_DIRS}/boost/math/quadrature/naive_monte_carlo.hpp)
if(boost_int_include)
  list(APPEND sources IntegratorNaive.cpp)
endif()

#----- build the object

cepgen_build(CepGenBoost SOURCES ${sources}
    EXT_LIBS ${Boost_LIBRARIES}
    EXT_HEADERS ${Boost_INCLUDE_DIRS}
    INSTALL_COMPONENT boost)
