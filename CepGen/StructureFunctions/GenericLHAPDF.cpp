#include "GenericLHAPDF.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  namespace SF
  {
    constexpr std::array<short,6> GenericLHAPDF::qtimes3_;

    GenericLHAPDF::Parameterisation::Parameterisation() :
      num_flavours( 4 ), pdf_set( "cteq6" ), pdf_member( 0 )
    {}

    GenericLHAPDF::Parameterisation
    GenericLHAPDF::Parameterisation::cteq6()
    {
      Parameterisation p;
      p.num_flavours = 4;
      p.pdf_set = "cteq6";
      p.pdf_member = 0;
      return p;
    }

    GenericLHAPDF::GenericLHAPDF( const Parameterisation& param ) :
      StructureFunctions( Type::GenericLHAPDF ), params( param ), initialised_( false )
    {}

    GenericLHAPDF::GenericLHAPDF( const char* set ) :
      StructureFunctions( Type::GenericLHAPDF ), initialised_( false )
    {
      params.pdf_set = set;
      initialise();
    }

    void
    GenericLHAPDF::initialise()
    {
      if ( initialised_ )
        return;
#ifdef LIBLHAPDF
      std::string lhapdf_version;
#  if LHAPDF_MAJOR_VERSION == 6
      try {
        pdf_set_ = LHAPDF::PDFSet( params.pdf_set );
        pdf_set_.mkPDFs<std::unique_ptr<LHAPDF::PDF> >( pdfs_ );
      } catch ( const LHAPDF::Exception& e ) {
        throw CG_FATAL( "GenericLHAPDF" )
          << "Caught LHAPDF exception:\n\t"
          << e.what();
      }
      lhapdf_version = LHAPDF::version();
#  else
      LHAPDF::initPDFSet( params.pdf_set, LHAPDF::LHGRID, params.pdf_member );
      lhapdf_version = LHAPDF::getVersion();
#  endif
      CG_INFO( "GenericLHAPDF" ) << "LHAPDF structure functions evaluator successfully built.\n"
        << " *) LHAPDF version: " << lhapdf_version << "\n"
        << " *) number of flavours: " << params.num_flavours << "\n"
        << " *) PDF set: " << params.pdf_set << "\n"
        << " *) PDF member: " << params.pdf_member;
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
        const double prefactor = 1./9.*qtimes3_[i]*qtimes3_[i];
#  if LHAPDF_MAJOR_VERSION == 6
        const double xq = pdfs_[params.pdf_member]->xfxQ2( i, xbj, q2 );
        const double xqbar = pdfs_[params.pdf_member]->xfxQ2( -i, xbj, q2 );
#  else
        const double xq = LHAPDF::xfx( xbj, q2, i+1 );
        const double xqbar = LHAPDF::xfx( xbj, q2, -( i+1 ) );
#  endif
        F2 += prefactor*( xq+xqbar );
      }
#else
      throw CG_FATAL( "GenericLHAPDF" ) << "LHAPDF is not liked to this instance!";
#endif

      return *this;
    }
  }
}
