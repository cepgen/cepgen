file(GLOB cards_sources Handler.cpp LpairHandler.cpp CommandLineHandler.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

cepgen_build(CepGenCards SOURCES ${cards_sources}
    INSTALL_COMPONENT lib)

