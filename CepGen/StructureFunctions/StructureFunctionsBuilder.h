#ifndef CepGen_StructureFunctions_StructureFunctionsBuilder_h
#define CepGen_StructureFunctions_StructureFunctionsBuilder_h

#include <memory>

namespace CepGen
{
  class StructureFunctions;
  namespace SF { enum class Type; }
  /// Helper class to generate any supported set of structure functions
  class StructureFunctionsBuilder
  {
    public:
      StructureFunctionsBuilder() {}
      ~StructureFunctionsBuilder() {}

      /// Build structure functions from the modelling type
      static std::shared_ptr<StructureFunctions> get( const SF::Type& );
      /// Build structure functions from the modelling name
      static std::shared_ptr<StructureFunctions> get( const char* );
  };
}

#endif
