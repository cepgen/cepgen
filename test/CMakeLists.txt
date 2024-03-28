cepgen_test_category(NAME "General"
                     SOURCES *.cc
                     FSOURCES *.f)

add_subdirectory(benchmarks)
add_subdirectory(generator)
add_subdirectory(physics)
add_subdirectory(utils)
add_subdirectory(validation)

set(CEPGEN_UNSTABLE_TESTS ${CEPGEN_UNSTABLE_TESTS} PARENT_SCOPE)
