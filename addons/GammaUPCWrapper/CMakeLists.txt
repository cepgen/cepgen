#--- searching for gammaUPC
set(GAMMAUPC_DIRS $ENV{GAMMAUPC_PATH} $ENV{GAMMAUPC_DIR} ${GAMMAUPC_DIR})
find_path(gammaUPC NAMES ElasticPhotonPhotonFlux.f90 HINTS ${GAMMAUPC_DIRS} PATH_SUFFIXES src)
if (NOT gammaUPC)
    return()
endif ()

file(GLOB gammaUPC_SOURCES ${gammaUPC}/*.f90)
message(STATUS "gammaUPC ${gammaUPC} sources: ${gammaUPC_SOURCES}")

#----- build the objects
cepgen_build(CepGenGammaUPC
        SOURCES src/*.cpp ${gammaUPC_SOURCES}
        OPTIONS "-Wa,-Wno-tabs -Wno-unused-variable -Wno-unused-function -Wno-maybe-uninitialized"
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
        COMPONENT gammaUPC)
cpack_add_component(gammaUPC
        DISPLAY_NAME " CepGen gammaUPC wrappers library"
        DESCRIPTION "Collection of CepGen wrappers to the gammaUPC library"
        DEPENDS lib)
