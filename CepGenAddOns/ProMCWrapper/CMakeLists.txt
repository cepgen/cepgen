#--- linking with ProMC

set(PROMC_DIRS $ENV{PROMC_DIR} $ENV{PROMC} /usr /usr/local /usr/local/promc)
find_library(PROMC_LIB NAMES promc HINTS ${PROMC_DIRS} PATH_SUFFIXES lib)
find_path(PROMC_INCLUDE NAMES ProMCBook.h HINTS ${PROMC_DIRS} PATH_SUFFIXES src)
find_path(PROMC_EXT_INCLUDE NAMES CBook HINTS ${PROMC_DIRS} PATH_SUFFIXES include)

if(PROMC_LIB)
  message(STATUS "ProMC found in ${PROMC_LIB}")
else()
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources *.cpp)

#----- build the objects

add_library(CepGenProMC SHARED ${sources})
target_link_libraries(CepGenProMC PRIVATE ${PROMC_LIB})
target_include_directories(CepGenProMC PRIVATE ${PROMC_INCLUDE} ${PROMC_EXT_INCLUDE})

#----- installation rules

install(TARGETS CepGenProMC DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT promc)
