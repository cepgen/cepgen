cepgen_test_category(NAME "Generator"
                     SOURCES *.cc
                     FSOURCES *.f)
set_tests_properties(ktfactor_proc_hadronisation PROPERTIES DISABLED TRUE)  # requires an algorithm as argument

#--- one instance of the common kt-factorised hadronisation test per algorithm registered by the add-ons
get_property(hadronisers GLOBAL PROPERTY CEPGEN_HADRONISERS_TO_TEST)
foreach(algo ${hadronisers})
    message(STATUS "... added test ${algo}_kt_hadr")
    add_test(NAME ${algo}_kt_hadr
             COMMAND ktfactor_proc_hadronisation -H ${algo}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR})  # add-ons libraries, CepGenAddOns.txt, PDG table and README live there
    set_tests_properties(${algo}_kt_hadr PROPERTIES
                         ENVIRONMENT "CEPGEN_PATH=${PROJECT_SOURCE_DIR}"
                         LABELS "hadronisation")
endforeach()
