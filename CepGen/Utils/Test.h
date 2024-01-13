#ifndef CepGen_Utils_Test_h
#define CepGen_Utils_Test_h

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace test {
    bool debug = false;
    double failure_tolerance = 0.;
    const double base_precision = 1.e-3;
    double precision = base_precision;
    size_t num_total = 0;
    size_t num_passed = 0;
  }  // namespace test
}  // namespace cepgen

#define CG_FAILED(name)                                                                                            \
  CG_LOG << cepgen::utils::colourise("FAILED ", cepgen::utils::Colour::red, cepgen::utils::Modifier::bold) << name \
         << "!"
#define CG_PASSED(name) CG_LOG << cepgen::utils::colourise("Passed ", cepgen::utils::Colour::green) << name << "."

#define CG_TEST_DEBUG(debugging) cepgen::test::debug = debugging
#define CG_TEST_SET_FAILURE_TOLERANCE_RATE(tolerance) cepgen::test::failure_tolerance = tolerance
#define CG_TEST_SET_PRECISION(precis) cepgen::test::precision = precis
#define CG_TEST_RESET_PRECISION() cepgen::test::precision = cepgen::test::base_precision

#define CG_TEST(test_cond, name)                                                                                     \
  {                                                                                                                  \
    if (cepgen::test::debug)                                                                                         \
      CG_LOG << cepgen::utils::colourise("TEST INFO", cepgen::utils::Colour::magenta, cepgen::utils::Modifier::bold) \
             << " " << cepgen::utils::colourise(name, cepgen::utils::Colour::magenta) << "\n"                        \
             << "\tcondition: " << cepgen::utils::boldify(#test_cond) << ".";                                        \
    if (!(test_cond))                                                                                                \
      CG_FAILED(name);                                                                                               \
    else {                                                                                                           \
      CG_PASSED(name);                                                                                               \
      cepgen::test::num_passed++;                                                                                    \
    }                                                                                                                \
    cepgen::test::num_total++;                                                                                       \
  }

#define CG_TEST_EQUAL(var1, var2, name)                                                                              \
  {                                                                                                                  \
    if (cepgen::test::debug)                                                                                         \
      CG_LOG << cepgen::utils::colourise("TEST INFO", cepgen::utils::Colour::magenta, cepgen::utils::Modifier::bold) \
             << " " << cepgen::utils::colourise(name, cepgen::utils::Colour::magenta) << "\n"                        \
             << "\tvariable 1(" << cepgen::utils::boldify(#var1) << "): " << var1 << "\n"                            \
             << "\tvariable 2(" << cepgen::utils::boldify(#var2) << "): " << var2 << ".";                            \
    if ((var1) != (var2))                                                                                            \
      CG_FAILED(name) << " " << var1 << " != " << var2 << ".";                                                       \
    else {                                                                                                           \
      CG_PASSED(name);                                                                                               \
      cepgen::test::num_passed++;                                                                                    \
    }                                                                                                                \
    cepgen::test::num_total++;                                                                                       \
  }

#define CG_TEST_EQUIV(var1, var2, name)                                                                        \
  {                                                                                                            \
    if (std::fabs((var1) - (var2)) > cepgen::test::precision)                                                  \
      CG_FAILED(name) << " " << var1 << " is not within " << cepgen::test::precision << " of " << var2 << "."; \
    else {                                                                                                     \
      CG_PASSED(name);                                                                                         \
      cepgen::test::num_passed++;                                                                              \
    }                                                                                                          \
    cepgen::test::num_total++;                                                                                 \
  }

#define CG_TEST_UNCERT(diff, unc, nsigma, name)                                                                      \
  {                                                                                                                  \
    if (cepgen::test::debug)                                                                                         \
      CG_LOG << cepgen::utils::colourise("TEST INFO", cepgen::utils::Colour::magenta, cepgen::utils::Modifier::bold) \
             << " " << cepgen::utils::colourise(name, cepgen::utils::Colour::magenta) << "\n"                        \
             << "\tdifference: " << diff << ", sigma: " << unc << " = " << diff / unc << " * sigma "                 \
             << ((diff > nsigma * unc) ? ">" : "<") << " " << nsigma << " * sigma.";                                 \
    if (diff > nsigma * unc)                                                                                         \
      CG_FAILED(name) << " difference " << diff << " is not within " << nsigma << " sigmas=" << unc << ".";          \
    else {                                                                                                           \
      CG_PASSED(name);                                                                                               \
      cepgen::test::num_passed++;                                                                                    \
    }                                                                                                                \
    cepgen::test::num_total++;                                                                                       \
  }

#define CG_TEST_EXCEPT(sequence, name)                                                                                \
  if (cepgen::test::debug)                                                                                            \
    CG_LOG << cepgen::utils::colourise("TEST INFO", cepgen::utils::Colour::magenta, cepgen::utils::Modifier::bold)    \
           << " " << cepgen::utils::colourise(name, cepgen::utils::Colour::magenta) << "\n"                           \
           << "\tsequence: " << cepgen::utils::boldify(#sequence) << ".";                                             \
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
      cepgen::test::num_passed++;                                                                                     \
    }                                                                                                                 \
    cepgen::test::num_total++;                                                                                        \
  }

#define CG_TEST_SUMMARY                                                                                            \
  {                                                                                                                \
    {                                                                                                              \
      auto col = cepgen::utils::Colour::yellow;                                                                    \
      if (cepgen::test::num_passed == cepgen::test::num_total)                                                     \
        col = cepgen::utils::Colour::green;                                                                        \
      else if ((cepgen::test::failure_tolerance > 0. &&                                                            \
                cepgen::test::num_total - cepgen::test::num_passed >                                               \
                    cepgen::test::failure_tolerance * cepgen::test::num_total) ||                                  \
               cepgen::test::num_passed < 0.1 * cepgen::test::num_total)                                           \
        col = cepgen::utils::Colour::red;                                                                          \
      CG_LOG << cepgen::utils::colourise(std::to_string(cepgen::test::num_passed) + " out of " +                   \
                                             cepgen::utils::s("test", cepgen::test::num_total, true) + " passed.", \
                                         col);                                                                     \
    }                                                                                                              \
    return (cepgen::test::num_total - cepgen::test::num_passed <=                                                  \
            cepgen::test::failure_tolerance * cepgen::test::num_total)                                             \
               ? 0                                                                                                 \
               : cepgen::test::num_total - cepgen::test::num_passed;                                               \
  }

#endif
