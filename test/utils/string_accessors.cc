#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  bool verbose;
  cepgen::ArgumentsParser(argc, argv).addOptionalArgument("verbose", "verbose mode", &verbose, false).parse();
  CG_TEST_DEBUG(verbose);
  {
    auto str = "This {is} a text with {sub-strings}";
    const auto substr = cepgen::utils::between(str, "{", "}");
    CG_TEST_EQUAL(substr.size(), 2, "number of substrings");
    const auto exp_substr = vector<string>{"is", "sub-strings"};
    CG_TEST_EQUAL(substr, exp_substr, "substrings content");
  }
  CG_TEST_SUMMARY;
}
