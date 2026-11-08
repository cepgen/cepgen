include(FetchContent)
FetchContent_Declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench.git
    GIT_TAG v4.3.11
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(nanobench)

cepgen_test_category(NAME "Benchmarks"
                     BSOURCES *.cc
                     LIBRARIES nanobench
                     PREPEND benchmark_
                     FSOURCES *.f)
