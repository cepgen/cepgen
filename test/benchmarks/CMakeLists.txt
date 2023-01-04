find_package(benchmark QUIET)
if(NOT benchmark_FOUND)
  return()
endif()

cepgen_test_category(NAME "Benchmarks"
                     BSOURCES *.cc
                     EXT_LIBS benchmark::benchmark
                     PREPEND benchmark_
                     FSOURCES *.f)
