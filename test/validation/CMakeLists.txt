file(MAKE_DIRECTORY ${CEPGEN_PATH}/validation)

cepgen_test_category(NAME "Validation"
    SOURCES *.cc
    WORKING_DIRECTORY ${CEPGEN_PATH}
    PREPEND validation_)
