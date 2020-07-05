#ifndef CepGen_Module_NamedModule_h
#define CepGen_Module_NamedModule_h

#include "CepGen/Core/ParametersList.h"

namespace cepgen
{
  template<typename T=std::string>
  class NamedModule
  {
    public:
      explicit NamedModule( const ParametersList& params ) :
        params_( params ),
        name_( params.name<T>( "<invalid>" ) ),
        description_( params.get<std::string>( "description" ) ) {}

      /// Module unique name
      const T& name() const { return name_; }
      /// Module description
      const std::string& description() const { return description_; }
      /// Module user-defined parameters
      inline const ParametersList& parameters() const { return params_; }

    protected:
      /// Set of parameters to steer this output module
      const ParametersList params_;
      /// Module unique name
      const T name_;
      /// Module description
      const std::string description_;
  };
}

#endif

