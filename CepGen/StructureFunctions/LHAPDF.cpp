#include "CepGen/StructureFunctions/LHAPDF.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

namespace CepGen
{
  namespace SF
  {
    constexpr std::array<short,6> LHAPDF::qtimes3_;

    LHAPDF::Parameterisation::Parameterisation() :
      num_flavours( 4 ), pdf_set( "cteq6" ), pdf_member( 0 )
    {}

    LHAPDF::Parameterisation
    LHAPDF::Parameterisation::cteq6()
    {
      Parameterisation p;
      p.num_flavours = 4;
      p.pdf_set = "cteq6";
      p.pdf_member = 0;
      return p;
    }

    LHAPDF::LHAPDF( const Parameterisation& param ) :
      StructureFunctions( Type::LHAPDF ), params( param ), initialised_( false )
    {}

    LHAPDF::LHAPDF( const char* set ) :
      StructureFunctions( Type::LHAPDF ), initialised_( false )
    {
      params.pdf_set = set;
      initialise();
    }

    void
    LHAPDF::initialise()
    {
      if ( initialised_ )
        return;
#ifdef LIBLHAPDF
      std::string lhapdf_version, pdf_description, pdf_type;
#  if LHAPDF_MAJOR_VERSION == 6
      try {
        pdf_set_ = ::LHAPDF::PDFSet( params.pdf_set );
        pdf_set_.mkPDFs<std::unique_ptr<::LHAPDF::PDF> >( pdfs_ );
        lhapdf_version = ::LHAPDF::version();
        pdf_description = pdf_set_.description();
        pdf_type = pdfs_[params.pdf_member]->type();
      } catch ( const ::LHAPDF::Exception& e ) {
        throw CG_FATAL( "LHAPDF" )
          << "Caught LHAPDF exception:\n\t"
          << e.what();
      }
#  else
      ::LHAPDF::initPDFSet( params.pdf_set, ::LHAPDF::LHGRID, params.pdf_member );
      lhapdf_version = ::LHAPDF::getVersion();
      pdf_description = ::LHAPDF::getDescription();
#  endif
      replace_all( pdf_description, ". ", ".\n  " );
      CG_INFO( "LHAPDF" ) << "LHAPDF structure functions evaluator successfully built.\n"
        << " *) LHAPDF version: " << lhapdf_version << "\n"
        << " *) number of flavours: " << params.num_flavours << "\n"
        << " *) PDF set: " << params.pdf_set << "\n"
        << ( pdf_description.empty() ? "" : "  "+pdf_description+"\n" )
        << " *) PDF member: " << params.pdf_member << ( pdf_type.empty() ? "" : " ("+pdf_type+")" );
      initialised_ = true;
#else
      throw CG_FATAL( "LHAPDF" ) << "LHAPDF is not liked to this instance!";
#endif
    }

    LHAPDF&
    LHAPDF::operator()( double q2, double xbj )
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
      throw CG_FATAL( "LHAPDF" ) << "LHAPDF is not liked to this instance!";
#endif

      return *this;
    }
  }
}
