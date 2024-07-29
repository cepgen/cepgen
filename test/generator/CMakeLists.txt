cepgen_test_category(NAME "Generator"
                     SOURCES *.cc
                     FSOURCES *.f)
set_tests_properties(ktfactor_proc_hadronisation PROPERTIES DISABLED TRUE)  # requires an algorithm as argument
