if(NOT ROOT_FOUND)
  message(FATAL_ERROR "Trying to build ROOT wrapper library while no ROOT installation was found on system!")
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources *.cpp)

#----- build the objects

add_library(CepGenRoot SHARED ${sources})
target_compile_definitions(CepGenRoot PUBLIC ${ROOT_DEFINITIONS})
target_link_libraries(CepGenRoot PUBLIC ${ROOT_LIBRARIES})
target_include_directories(CepGenRoot PUBLIC ${ROOT_INCLUDE_DIRS})

#----- installation rules

install(TARGETS CepGenRoot DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT root)
