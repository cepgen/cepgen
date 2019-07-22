#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Hadronisers/HadronisersHandler.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Core/ParametersList.h"

namespace cepgen
{
  namespace hadr
  {
    /// Interface to the Tauola decay routine
    class TauolaFilter : public GenericHadroniser
    {
      public:
        explicit TauolaFilter( const ParametersList& );

        void setParameters( const Parameters& ) override {}
        inline void readString( const char* param ) override;
        void init() override {}
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override {}

      private:

    };

    TauolaFilter::TauolaFilter( const ParametersList& plist ) :
      GenericHadroniser( plist, "tauola" )
    {}

    void
    TauolaFilter::readString( const char* param )
    {
    }

    bool
    TauolaFilter::run( Event& ev, double& weight, bool full )
    {
      weight = 1.;

      return true;
    }
  }
}

// register hadroniser and define alias
REGISTER_HADRONISER( tauola, TauolaFilter )

