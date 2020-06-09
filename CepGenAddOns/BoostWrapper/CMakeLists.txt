if(NOT Boost_FOUND)
  message(FATAL_ERROR "Trying to build Boost wrapper library while no Boost installation was found on system!")
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB cards_sources BoostTreeHandler.cpp)

#----- build the objects

add_library(CepGenBoost SHARED ${cards_sources})
target_link_libraries(CepGenBoost PRIVATE ${Boost_LIBRARIES})
target_include_directories(CepGenBoost PRIVATE ${Boost_INCLUDE_DIRS})
list(APPEND CEPGEN_ADDONS_EXT CepGenBoost)
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenBoost DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT boost)
