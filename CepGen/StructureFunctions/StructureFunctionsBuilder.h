#ifndef CepGen_StructureFunctions_StructureFunctionsBuilder_h
#define CepGen_StructureFunctions_StructureFunctionsBuilder_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

namespace CepGen
{
  /// Helper class to generate any supported set of structure functions
  class StructureFunctionsBuilder
  {
    public:
      StructureFunctionsBuilder() {}
      ~StructureFunctionsBuilder() {}

      /// Build structure functions from the modelling type
      static StructureFunctions get( const StructureFunctions::Type&, double q2, double xbj );
      /// Build structure functions from the modelling name
      static StructureFunctions get( const char*, double q2, double xbj );
  };
}

#endif
