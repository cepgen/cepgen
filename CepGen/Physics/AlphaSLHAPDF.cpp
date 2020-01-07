#include "CepGen/Physics/AlphaS.h"

#include "LHAPDF/LHAPDF.h"

#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
# define LHAPDF_GE_6 1
#endif

namespace cepgen
{
  class AlphaSLHAPDF : public AlphaS
  {
    public:
      explicit AlphaSLHAPDF( const ParametersList& params )
#ifdef LHAPDF_GE_6
      : lhapdf_( LHAPDF::mkPDF(
        params.get<std::string>( "pdfSet", "cteq6" ),
        params.get<int>( "pdfMember", 0 ) ) )
#endif
      {
#ifndef LHAPDF_GE_6
        LHAPDF::initPDFSet(
          params.get<std::string>( "pdfSet" ),
          LHAPDF::LHGRID,
          params.get<int>( "pdfMember", 0 ) );
#endif
      }

      double operator()( double q ) const override {
#ifdef LHAPDF_GE_6
        return lhapdf_->alphasQ( q );
#else
        return LHAPDF::alphasPDF( q );
#endif
      }

    private:
#ifdef LHAPDF_GE_6
      std::unique_ptr<LHAPDF::PDF> lhapdf_;
#endif
  };
}

REGISTER_ALPHAS_MODULE( "lhapdf", AlphaSLHAPDF )
