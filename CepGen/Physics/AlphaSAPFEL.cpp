#include "CepGen/Physics/AlphaS.h"

#include "APFEL/APFEL.h"

namespace cepgen
{
  class AlphaSAPFEL : public AlphaS
  {
    public:
      explicit AlphaSAPFEL( const ParametersList& params ) :
        order_( params.get<int>( "order", 2 ) ),
        q0_( params.get<double>( "q0", 1. ) ),
        qmax_( params.get<double>( "qmax", 100. ) ) {
        APFEL::SetPerturbativeOrder( order_ );
        APFEL::InitializeAPFEL();
        APFEL::EvolveAPFEL( q0_, qmax_ );
      }

      double operator()( double q ) const override {
        return APFEL::AlphaQCD( q );
      }

    private:
      int order_;
      double q0_, qmax_;
  };
}

REGISTER_ALPHAS_MODULE( "apfel", AlphaSAPFEL )
