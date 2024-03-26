#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  bool verbose;
  cepgen::ArgumentsParser(argc, argv).addOptionalArgument("verbose", "verbose mode", &verbose, false).parse();
  CG_TEST_DEBUG(verbose);
  {
    const auto kebab_str = "this-is-a-kebab-string";
    CG_TEST_EQUAL(cepgen::utils::toCamelCase(kebab_str), "thisIsAKebabString", "kebab -> camel");
  }
  {
    const auto snake_str = "this_is_a_snake_case_string";
    CG_TEST_EQUAL(cepgen::utils::toCamelCase(snake_str), "thisIsASnakeCaseString", "snake -> camel");
  }
  {
    const auto screaming_snake_str = "THIS_IS_A_SCREAMING_SNAKE_CASE_STRING";
    CG_TEST_EQUAL(
        cepgen::utils::toCamelCase(screaming_snake_str), "thisIsAScreamingSnakeCaseString", "screaming snake -> camel");
  }
  {
    const auto camel_str = "thisIsACamelCaseString";
    CG_TEST_EQUAL(cepgen::utils::toCamelCase(camel_str), camel_str, "camel -> camel");
  }
  CG_TEST_SUMMARY;
}
