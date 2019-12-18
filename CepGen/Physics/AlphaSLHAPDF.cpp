#include "CepGen/Physics/AlphaS.h"

#include "LHAPDF/LHAPDF.h"

namespace cepgen
{
  class AlphaSLHAPDF : public AlphaS
  {
    public:
      explicit AlphaSLHAPDF( const ParametersList& params ) :
      lhapdf_( LHAPDF::mkPDF(
        params.get<std::string>( "pdfSet", "cteq6" ),
        params.get<int>( "pdfMember", 0 ) ) ) {
      }

      double operator()( double q ) const override {
        return lhapdf_->alphasQ( q );
      }

    private:
      std::unique_ptr<LHAPDF::PDF> lhapdf_;
  };
}

REGISTER_ALPHAS_MODULE( "lhapdf", AlphaSLHAPDF )
