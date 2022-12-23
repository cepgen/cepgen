cepgen_test_category(NAME "General"
                     SOURCES *.cc
                     FSOURCES *.f)

add_subdirectory(generator)
add_subdirectory(physics)
add_subdirectory(utils)
