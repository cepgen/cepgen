#include "GenericLHAPDF.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  namespace SF
  {
    constexpr std::array<double,6> GenericLHAPDF::qtimes3_;

    GenericLHAPDF::Parameterisation
    GenericLHAPDF::Parameterisation::cteq6()
    {
      Parameterisation p;
      p.num_flavours = 4;
      p.pdf_set = "cteq6";
      return p;
    }

    GenericLHAPDF::GenericLHAPDF( const Parameterisation& param ) :
      StructureFunctions( Type::GenericLHAPDF ), params( param ), initialised_( false )
    {}

    GenericLHAPDF::GenericLHAPDF( const char* set ) :
      StructureFunctions( Type::GenericLHAPDF )
    {
      params.pdf_set = set;
      initialise();
    }

    void
    GenericLHAPDF::initialise()
    {
#ifdef LIBLHAPDF
#  if LHAPDF_MAJOR_VERSION == 6
      pdf_set_ = LHAPDF::PDFSet( params.pdf_set );
      pdfs_ = pdf_set_.mkPDFs();
#  else
      LHAPDF::initPDFSet( params.pdf_set, LHAPDF::LHGRID, 0 );
#  endif
      CG_INFO( "GenericLHAPDF" ) << "LHAPDF structure functions evaluator successfully built.\n"
        << " *) number of flavours: " << params.num_flavours << "\n"
        << " *) PDFset: " << params.pdf_set;
      initialised_ = true;
#else
      throw CG_FATAL( "GenericLHAPDF" ) << "LHAPDF is not liked to this instance!";
#endif
    }

    GenericLHAPDF&
    GenericLHAPDF::operator()( double q2, double xbj )
    {
      std::pair<double,double> nv = { q2, xbj };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      F2 = 0.;
      if ( params.num_flavours == 0 || params.num_flavours > 6 )
        return *this;

      if ( !initialised_ )
        initialise();
#ifdef LIBLHAPDF
      for ( int i = 0; i < params.num_flavours; ++i ) {
#  if LHAPDF_MAJOR_VERSION == 6
        const double xq = pdfs_[0]->xfxQ2( i, xbj, q2 );
        const double xqbar = pdfs_[0]->xfxQ2( -i, xbj, q2 );
#  else
        const double xq = LHAPDF::xfx( xbj, q2, i+1 );
        const double xqbar = LHAPDF::xfx( xbj, q2, -i-1 );
#  endif
        F2 += qtimes3_[i]*qtimes3_[i]/9. * ( xq + xqbar );
      }
#else
      throw CG_FATAL( "GenericLHAPDF" ) << "LHAPDF is not liked to this instance!";
#endif

      return *this;
    }
  }
}
