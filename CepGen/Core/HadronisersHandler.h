#ifndef CepGen_Processes_HadronisersHandler_h
#define CepGen_Processes_HadronisersHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#define BUILDERNM( obj ) obj ## Builder
#define STRINGIFY( name ) #name
#define REGISTER_HADRONISER( name, obj ) \
  struct BUILDERNM( name ) { \
    BUILDERNM( name )() { cepgen::HadronisersHandler::get().registerModule( STRINGIFY( name ), new obj ); } }; \
  static BUILDERNM( name ) g ## name;

namespace cepgen
{
  namespace hadr
  {
    class HadronisersHandler : public ModuleFactory<GenericHadroniser>
    {
      public:
        std::unique_ptr<GenericHadroniser> build( const std::string& name, const ParametersList& ) const override {
          if ( map_.count( name ) == 0 )
            throw std::logic_error( "Failed to retrieve a hadroniser with name \""+name+"\"!" );
          return std::unique_ptr<GenericHadroniser>( new GenericHadroniser( *map_.at( name ) ) );
        }
    };
  }
}

#endif

