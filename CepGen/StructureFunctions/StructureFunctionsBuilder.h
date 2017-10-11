#ifndef CepGen_StructureFunctions_StructureFunctionsBuilder_h
#define CepGen_StructureFunctions_StructureFunctionsBuilder_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

namespace CepGen
{
  class StructureFunctionsBuilder
  {
    public:
      StructureFunctionsBuilder() {}
      ~StructureFunctionsBuilder() {}

      static StructureFunctions get( const StructureFunctions::Type&, double q2, double xbj );
      static StructureFunctions get( const char*, double q2, double xbj );
  };
}

#endif
