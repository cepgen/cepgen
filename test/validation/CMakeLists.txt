file(MAKE_DIRECTORY $ENV{CEPGEN_PATH}/validation)

cepgen_test_category(NAME "Validation"
    SOURCES *.cc
    PREPEND validation_)
