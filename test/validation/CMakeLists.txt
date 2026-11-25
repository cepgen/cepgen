file(MAKE_DIRECTORY ${CEPGEN_PATH}/validation)

cepgen_test_category(NAME "Validation"
    SOURCES *.cc
    WORKING_DIRECTORY ${CEPGEN_PATH}
    PREPEND validation_)
add_test(NAME validation_compare_dilepton_kt   COMMAND validation_compare_dilepton --processes "lpair;pptoff<method:0;pptoff<method:1;pptoff<method:2" --sqrts 14000 --filename validation/comparison_dilepton_kt_ --kt)
add_test(NAME validation_compare_dilepton_coll COMMAND validation_compare_dilepton --processes "lpair;pptoff<method:0;mg5_aMC" --sqrts 14000 --filename validation/comparison_dilepton_coll_)
set_tests_properties(validation_compare_dilepton PROPERTIES DISABLED TRUE)
