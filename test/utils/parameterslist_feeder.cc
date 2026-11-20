#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  {
    const string fed = "test/of/key:value";
    cepgen::ParametersList plist;
    plist.feed(fed);
    CG_TEST_EQUAL(plist.get<cepgen::ParametersList>("test").get<cepgen::ParametersList>("of").get<std::string>("key"),
                  "value",
                  "parameters list chain");
    CG_DEBUG("main") << "Resulting parameters list: " << plist << ".";
    CG_TEST_EQUAL(plist.serialise(), fed, "parameters list serialisation");
  }
  {
    cepgen::ParametersList plist;
    plist.feed("foo:3.14").feed("bar:2").feed("baz:2e3");
    CG_DEBUG("main") << "Resulting parameters list: " << plist << ".";
    CG_TEST_EQUAL(plist.get<double>("foo"), 3.14, "float parsing");
    CG_TEST_EQUAL(plist.get<int>("bar"), 2, "integer parsing");
    CG_TEST_EQUAL(plist.get<double>("baz"), 2000., "float (from Ee notation");
    plist.feed("bat:5E10").feed("foo:42");
    CG_TEST_EQUAL(plist.get<double>("bat"), 5e10, "float (from re-parsing)");
    CG_TEST_EQUAL(plist.get<int>("foo"), 42, "integer (from re-parsing)");
    CG_TEST_EQUAL(plist.get<double>("foo", -1.), 42., "integer as float (from re-parsing)");
  }
  {
    auto fed = "this/is/a:test,this/works:true,that/{one:42,other:3.141592}";
    cepgen::ParametersList plist;
    plist.feed(fed);
    const auto re_fed = cepgen::ParametersList().feed(plist.serialise());
    CG_DEBUG("main") << "\n"
                     << "Fed string: " << fed << "\n"
                     << "Fed parameters list:\n"
                     << cepgen::ParametersDescription(plist) << "\n"
                     << "Re-serialised string: " << plist.serialise() << "\n"
                     << "Re-fed parameters list:\n"
                     << cepgen::ParametersDescription(re_fed) << ".\n"
                     << "Diff:\n"
                     << cepgen::ParametersDescription(plist.diff(re_fed)) << ".";
    CG_TEST_EQUAL(re_fed, plist, "serialised parameters list parsing");
  }
  CG_TEST_EXCEPT([]() { cepgen::ParametersList().feed("invalid/string/{{fed:true}"); }, "parsing of an invalid string");

  {  // a key is only held with one type
    cepgen::ParametersList plist;
    plist.set("key", std::vector<double>{0., 1.}).set("key", cepgen::Limits{2., 3.});
    CG_TEST(!plist.has<std::vector<double> >("key"), "single-type storage after set");
    CG_TEST_EQUAL(plist.typeOf("key"), "Limits", "stored type after set");
    CG_TEST_EQUAL(plist.keys().size(), 1ul, "single key after set");
    plist.operator[]<int>("key") = 4;
    CG_TEST(!plist.has<cepgen::Limits>("key"), "single-type storage after operator[]");
    plist.set("conv", std::vector<int>{1, 2});
    CG_TEST(plist.operator[]<cepgen::Limits>("conv") == cepgen::Limits(1., 2.), "in-place conversion by operator[]");
    CG_TEST(!plist.has<std::vector<int> >("conv"), "single-type storage after in-place conversion");
    CG_TEST_EQUAL(plist.typeOf("other"), "", "stored type of a missing key");
    cepgen::ParametersList other;
    other.set("sub", cepgen::ParametersList().set("a", 1));
    other += cepgen::ParametersList().set("sub", 2);
    CG_TEST_EQUAL(other.typeOf("sub"), "int", "single-type storage after concatenation");
  }
  {  // lossless scalar conversions
    const auto plist = cepgen::ParametersList()
                           .set("int", 42)
                           .set("negint", -3)
                           .set("dbl_int", 3.)
                           .set("dbl", 3.5)
                           .set("ulong", 12ull)
                           .set("ulong_big", (1ull << 53) + 1)
                           .set("bool", true)
                           .set("str", "42"s);
    CG_TEST_EQUAL(plist.get<double>("int"), 42., "double <- int");
    CG_TEST_EQUAL(plist.get<double>("ulong"), 12., "double <- ulong");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<double>("ulong_big"); }, "double <- non-representable ulong");
    CG_TEST_EQUAL(plist.get<int>("dbl_int"), 3, "int <- integral double");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<int>("dbl"); }, "int <- non-integral double");
    CG_TEST_EQUAL(plist.get<int>("bool"), 1, "int <- bool");
    CG_TEST_EQUAL(plist.get<int>("ulong"), 12, "int <- ulong");
    CG_TEST_EQUAL(plist.get<unsigned long long>("int"), 42ull, "ulong <- int");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<unsigned long long>("negint"); }, "ulong <- negative int");
    CG_TEST_EQUAL(cepgen::ParametersList().set("b", 1).get<bool>("b"), true, "bool <- int (1)");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<bool>("int"); }, "bool <- int (not 0/1)");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<int>("str"); }, "int <- string");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<std::string>("int"); }, "string <- int");
    CG_TEST_EQUAL(plist.get<double>("missing", -1.), -1., "default value for a missing key");
    int cast_value = -1;
    CG_TEST(!plist.cast("dbl", cast_value), "non-throwing cast of non-convertible value");
    CG_TEST_EQUAL(cast_value, -1, "non-convertible value leaves variable untouched");
    CG_TEST(plist.cast("dbl_int", cast_value), "non-throwing cast of convertible value");
    CG_TEST_EQUAL(cast_value, 3, "convertible value cast");
    int value = -1;
    plist.fill("dbl_int", value);
    CG_TEST_EQUAL(value, 3, "fill with conversion");
    CG_TEST_EXCEPT(
        [&plist]() {
          int value;
          plist.fill("dbl", value);
        },
        "fill with a non-convertible value");
  }
  {  // limits and vectors conversions
    const auto plist = cepgen::ParametersList()
                           .set("vint", std::vector<int>{0, 1})
                           .set("vdbl", std::vector<double>{0.5, 2.})
                           .set("vdbl1", std::vector<double>{0.5})
                           .set("vdbl3", std::vector<double>{0., 1., 2.})
                           .set("lim", cepgen::Limits{1., 10.})
                           .set("empty", std::vector<int>{})
                           .set("vvdbl", std::vector<std::vector<double> >{{0., 1.}, {2., 3.}})
                           .set("vlim", std::vector<cepgen::Limits>{{0., 1.}});
    CG_TEST(plist.get<cepgen::Limits>("vint") == cepgen::Limits(0., 1.), "Limits <- vint (Python-like (0, 1))");
    CG_TEST(plist.get<cepgen::Limits>("vdbl") == cepgen::Limits(0.5, 2.), "Limits <- vfloat");
    CG_TEST(plist.get<cepgen::Limits>("vdbl1").min() == 0.5 && !plist.get<cepgen::Limits>("vdbl1").hasMax(),
            "Limits <- 1-element vfloat");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<cepgen::Limits>("vdbl3"); }, "Limits <- 3-element vfloat");
    cepgen::Limits lim;
    plist.fill("vint", lim);
    CG_TEST(lim == cepgen::Limits(0., 1.), "fill Limits <- vint");
    CG_TEST((plist.get<std::vector<double> >("vint") == std::vector<double>{0., 1.}), "vfloat <- vint");
    CG_TEST((plist.get<std::vector<double> >("lim") == std::vector<double>{1., 10.}), "vfloat <- Limits");
    CG_TEST((plist.get<std::vector<int> >("vdbl3") == std::vector<int>{0, 1, 2}), "vint <- integral vfloat");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<std::vector<int> >("vdbl"); }, "vint <- non-integral vfloat");
    CG_TEST(plist.get<std::vector<cepgen::ParametersList> >("empty").empty(), "any vector <- empty vector");
    CG_TEST(plist.get<std::vector<cepgen::Limits> >("vvdbl").at(1) == cepgen::Limits(2., 3.), "VLimits <- vvfloat");
    CG_TEST((plist.get<std::vector<std::vector<double> > >("vlim").at(0) == std::vector<double>{0., 1.}),
            "vvfloat <- VLimits");
    CG_TEST_EXCEPT([&plist]() { (void)plist.get<cepgen::ParametersList>("vint"); }, "Params <- vint");
    CG_TEST((plist.get<std::vector<cepgen::Limits> >("vdbl") == std::vector<cepgen::Limits>{{0.5, 2.}}),
            "VLimits <- single range");
    CG_TEST(cepgen::ParametersList().set("vl", std::vector<cepgen::Limits>{{1., 2.}}).get<cepgen::Limits>("vl") ==
                cepgen::Limits(1., 2.),
            "Limits <- 1-element VLimits");
    // legacy min/max attributes
    const auto legacy = cepgen::ParametersList().set("ptmin", 1.).set("ptmax", 5.);
    CG_TEST(legacy.get<cepgen::Limits>("pt") == cepgen::Limits(1., 5.), "Limits <- min/max attributes");
    auto merged = cepgen::ParametersList().set("pt", std::vector<double>{0., 10.});
    merged += cepgen::ParametersList().set("ptmax", 5.);
    CG_TEST(merged.get<cepgen::Limits>("pt") == cepgen::Limits(0., 5.), "concatenation of vfloat limits and max");
  }

  CG_TEST_SUMMARY;
}
