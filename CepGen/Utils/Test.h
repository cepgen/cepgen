#ifndef CepGen_Utils_Test_h
#define CepGen_Utils_Test_h

#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

size_t num_tests = 0;
size_t num_tests_passed = 0;

#define CG_FAILED(name)                                                                                            \
  CG_LOG << cepgen::utils::colourise("FAILED ", cepgen::utils::Colour::red, cepgen::utils::Modifier::bold) << name \
         << "!"
#define CG_PASSED(name) CG_LOG << cepgen::utils::colourise("Passed ", cepgen::utils::Colour::green) << name << "."

#define CG_TEST(test_cond, name) \
  {                              \
    if (!(test_cond))            \
      CG_FAILED(name);           \
    else {                       \
      CG_PASSED(name);           \
      num_tests_passed++;        \
    }                            \
    num_tests++;                 \
  }

#define CG_TEST_EQUAL(var1, var2, name)                        \
  {                                                            \
    if ((var1) != (var2))                                      \
      CG_FAILED(name) << " " << var1 << " != " << var2 << "."; \
    else {                                                     \
      CG_PASSED(name);                                         \
      num_tests_passed++;                                      \
    }                                                          \
    num_tests++;                                               \
  }

#define CG_TEST_EXCEPT(sequence, name)                                                                                \
  try {                                                                                                               \
    sequence();                                                                                                       \
    throw cepgen::Exception("", "this_test");                                                                         \
  } catch (const cepgen::Exception& exc) {                                                                            \
    if (exc.from() == "this_test")                                                                                    \
      CG_FAILED(name);                                                                                                \
    else {                                                                                                            \
      CG_PASSED(name) << " Resulting exception:\n"                                                                    \
                      << cepgen::utils::colourise(exc.message(),                                                      \
                                                  cepgen::utils::Colour::none,                                        \
                                                  cepgen::utils::Modifier::dimmed | cepgen::utils::Modifier::italic); \
      num_tests_passed++;                                                                                             \
    }                                                                                                                 \
    num_tests++;                                                                                                      \
  }

#define CG_TEST_SUMMARY                                                                                           \
  {                                                                                                               \
    {                                                                                                             \
      auto col = cepgen::utils::Colour::yellow;                                                                   \
      if (num_tests_passed == num_tests)                                                                          \
        col = cepgen::utils::Colour::green;                                                                       \
      else if (num_tests_passed < 0.1 * num_tests)                                                                \
        col = cepgen::utils::Colour::red;                                                                         \
      CG_LOG << cepgen::utils::colourise(                                                                         \
          std::to_string(num_tests_passed) + " out of " + cepgen::utils::s("test", num_tests, true) + " passed.", \
          col);                                                                                                   \
    }                                                                                                             \
    return (num_tests - num_tests_passed);                                                                        \
  }

#endif
