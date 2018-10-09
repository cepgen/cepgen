#include "CepGen/Core/Exception.h"

using namespace cepgen;

int main()
{
  throw Exception( "Test", Exception::Type::warning ) << "Haha";
  return 0;
}
