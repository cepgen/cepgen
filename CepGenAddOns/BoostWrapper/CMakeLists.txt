if(NOT Boost_FOUND)
  message(FATAL_ERROR "Trying to build Boost wrapper library while no Boost installation was found on system!")
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB cards_sources BoostTreeHandler.cpp)

#----- build the objects

add_library(CepGenBoost SHARED ${cards_sources})
target_link_libraries(CepGenBoost PUBLIC ${Boost_LIBRARIES})
target_include_directories(CepGenBoost PUBLIC ${Boost_INCLUDE_DIRS})

#----- installation rules

install(TARGETS CepGenBoost DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT boost)
