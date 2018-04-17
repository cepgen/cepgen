#include "CepGen/Core/Exception.h"

using namespace CepGen;

int main()
{
  throw Exception( "Test", kJustWarning ) << "Haha";
  return 0;
}
