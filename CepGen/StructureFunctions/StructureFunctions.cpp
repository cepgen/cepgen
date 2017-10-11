#include "StructureFunctions.h"
#include <iostream>

namespace CepGen
{
  /// Human-readable format of a structure function object
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctions& sf )
  {
    return os << "F2 = " << sf.F2 << ", FL = " << sf.FL;
  }
}
