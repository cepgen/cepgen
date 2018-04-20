#include "CepGen/Core/Exception.h"

using namespace CepGen;

int main()
{
  throw Exception( "Test", Exception::Type::warning ) << "Haha";
  return 0;
}
