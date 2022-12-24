cepgen_test_category(NAME "Utilities"
                     SOURCES *.cc
                     FSOURCES *.f)

add_test(NAME integrator_all COMMAND integrator -q -a)
