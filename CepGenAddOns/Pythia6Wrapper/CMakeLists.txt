set(PYTHIA6LIB)
if(PYTHIA6)
  message(STATUS "Pythia 6 found in ${PYTHIA6}")
  list(APPEND PYTHIA6LIB ${PYTHIA6})
  if(PYTHIA6DUMMY)
    message(STATUS "Pythia 6 dummy library found in ${PYTHIA6DUMMY}")
    list(APPEND PYTHIA6LIB ${PYTHIA6DUMMY})
  endif()
  message(STATUS "Pythia 8 found in ${PYTHIA8}")
else()
  file(GLOB PYTHIA6_SRC ${PROJECT_SOURCE_DIR}/External/pythia*.f)
  if(EXISTS $ENV{PYTHIA6_SRC})
    file(GLOB pythia6_src $ENV{PYTHIA6_SRC})
  elseif(EXISTS ${PYTHIA6_SRC})
    file(GLOB pythia6_src ${PYTHIA6_SRC})
  else()
    return()
  endif()
  message(STATUS "Pythia 6 source found in ${pythia6_src}")
  add_library(pythia6 SHARED ${pythia6_src})
  install(TARGETS pythia6 DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT pythia6)
  set_property(TARGET pythia6 PROPERTY Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-tabs -Wno-maybe-uninitialized -Wno-integer-division -Wno-unused-variable -Wno-unused-dummy-argument")
  list(APPEND PYTHIA6LIB pythia6)
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources *.cpp)

#----- build the objects

add_library(CepGenPythia6 SHARED ${sources})
target_link_libraries(CepGenPythia6 PRIVATE ${PYTHIA6LIB})
set_property(TARGET CepGenPythia6 PROPERTY CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
list(APPEND CEPGEN_ADDONS_EXT CepGenPythia6)
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenPythia6 DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT pythia6)

