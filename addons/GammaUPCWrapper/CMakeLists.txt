#--- searching for gammaUPC sources: either a standalone copy, or the one shipped with MadGraph5_aMC@NLO (>= 3.5.0)
set(GAMMAUPC_DIRS $ENV{GAMMAUPC_PATH} $ENV{gammaUPC_DIR} ${GAMMAUPC_DIR})
find_program(MADGRAPH_BIN NAMES mg5_aMC HINTS $ENV{MADGRAPH_DIR} /usr /usr/local PATH_SUFFIXES bin)
if (MADGRAPH_BIN)
    get_filename_component(mg_bin ${MADGRAPH_BIN} REALPATH)  # mg5_aMC may be a symbolic link
    get_filename_component(mg_dir ${mg_bin} DIRECTORY)
    get_filename_component(mg_dir ${mg_dir} DIRECTORY)
    list(APPEND GAMMAUPC_DIRS ${mg_dir}/Template/Common/Source/PDF/gammaUPC)
endif ()
find_path(gammaUPC NAMES ElasticPhotonPhotonFlux.f90 HINTS ${GAMMAUPC_DIRS} PATH_SUFFIXES src)
if (NOT gammaUPC)
    return()
endif ()

file(GLOB gammaUPC_SOURCES ${gammaUPC}/*.f90)
list(FILTER gammaUPC_SOURCES EXCLUDE REGEX ".*/test\\.f90$")  # standalone test program (with its own main)
message(STATUS "gammaUPC ${gammaUPC} sources: ${gammaUPC_SOURCES}")

#----- build the objects
cepgen_build(CepGenGammaUPC
        SOURCES src/*.cpp ${gammaUPC_SOURCES}
        TESTS test/*.cc
        OPTIONS "-Wa,-Wno-tabs -Wno-unused-variable -Wno-unused-function -Wno-maybe-uninitialized"
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
        COMPONENT gammaUPC)
cpack_add_component(gammaUPC
        DISPLAY_NAME " CepGen gammaUPC wrappers library"
        DESCRIPTION "Collection of CepGen wrappers to the gammaUPC library"
        DEPENDS lib)
