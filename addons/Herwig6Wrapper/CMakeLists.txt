#--- searching for Herwig 6
set(HERWIG6_DIRS ${CMAKE_CURRENT_BINARY_DIR} $ENV{HERWIG6_DIR} /usr /usr/local /opt/herwig6)
find_library(HERWIG6 NAMES herwig6 herwig HINTS ${HERWIG6_DIRS} PATH_SUFFIXES lib lib64)

set(HERWIG6LIB)
if (HERWIG6)
    list(APPEND HERWIG6LIB ${HERWIG6})
    find_library(HERWIG6DUMMY NAMES herwig_dummy herwig6_dummy HINTS ${HERWIG6_DIRS} PATH_SUFFIXES lib lib64)
    if (HERWIG6DUMMY)
        list(APPEND HERWIG6LIB ${HERWIG6DUMMY})
    endif ()
    find_library(HERWIG6PDFDUMMY NAMES herwig_pdfdummy herwig6_pdfdummy HINTS ${HERWIG6_DIRS} PATH_SUFFIXES lib lib64)
    if (HERWIG6PDFDUMMY)
        list(APPEND HERWIG6LIB ${HERWIG6PDFDUMMY})
    endif ()
else ()
    file(GLOB HERWIG6_SRC ${PROJECT_SOURCE_DIR}/External/herwig65*.f)
    file(GLOB HERWIG6_INC ${PROJECT_SOURCE_DIR}/External)
    if (EXISTS $ENV{HERWIG6_SRC})
        file(GLOB herwig6_src $ENV{HERWIG6_SRC})
    elseif (EXISTS ${HERWIG6_SRC})
        file(GLOB herwig6_src ${HERWIG6_SRC})
    else ()
        return()
    endif ()
    message(STATUS "Herwig 6 source found in ${herwig6_src}")
    add_library(herwig6 SHARED ${herwig6_src})
    target_include_directories(herwig6 PUBLIC ${HERWIG6_INC})
    target_compile_options(herwig6 PRIVATE -std=legacy)
    install(TARGETS herwig6 DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT herwig6)
    list(APPEND HERWIG6LIB herwig6)
endif ()

#----- build the object
cepgen_build(CepGenHerwig6
        SOURCES src/*.cpp
        TESTS test/*.cc
        LIBRARIES ${HERWIG6LIB}
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
        OPTIONS "-Wno-misleading-indentation"
        COMPONENT herwig6)
install(DIRECTORY CepGenHerwig6
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT herwig6)
