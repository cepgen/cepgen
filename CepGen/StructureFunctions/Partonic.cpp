#include "CepGen/StructureFunctions/Partonic.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

namespace cepgen
{
  namespace strfun
  {
    constexpr std::array<short,6> Partonic::Q_TIMES_3, Partonic::QUARK_PDGS;

    Partonic::Parameters::Parameters() :
      num_flavours( 4 ), pdf_set( "cteq6" ), pdf_code( 0l ), pdf_member( 0 ), mode( Mode::full )
    {}

    Partonic::Partonic( const Parameters& param ) :
      Parameterisation( Type::Partonic ), params( param ), initialised_( false )
    {}

    Partonic::Partonic( const char* set, unsigned short member, const Parameters::Mode& mode ) :
      Parameterisation( Type::Partonic ), initialised_( false )
    {
      params.pdf_set = set;
      params.pdf_member = member;
      params.mode = mode;
    }

    std::string
    Partonic::description() const
    {
      std::ostringstream os;
      os << "Partonic{" << params.pdf_set << ",m=" << params.pdf_member << ",mode=" << params.mode << "}";
      return os.str();
    }

    void
    Partonic::initialise()
    {
      if ( initialised_ )
        return;
#ifdef LIBLHAPDF
      std::string lhapdf_version, pdf_description, pdf_type;
#  if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
      try {
        //--- check if PDF code is set
        if ( params.pdf_code != 0l ) {
          auto pdf = LHAPDF::lookupPDF( params.pdf_code );
          if ( pdf.second != 0 )
            throw CG_FATAL( "Partonic" ) << "Failed to retrieve PDFset with id=" << params.pdf_code << "!";
          if ( !params.pdf_set.empty() && params.pdf_set != pdf.first )
            CG_WARNING( "Partonic" ) << "PDF set name changed from \"" << params.pdf_set << "\" to \"" << pdf.first << "\".";
          params.pdf_set = pdf.first;
        }
        pdf_set_ = LHAPDF::PDFSet( params.pdf_set );
        pdf_set_.mkPDFs<std::unique_ptr<LHAPDF::PDF> >( pdfs_ );
        lhapdf_version = LHAPDF::version();
        pdf_description = pdf_set_.description();
        pdf_type = pdfs_[params.pdf_member]->type();
      } catch ( const LHAPDF::Exception& e ) {
        throw CG_FATAL( "Partonic" )
          << "Caught LHAPDF exception:\n\t"
          << e.what();
      }
#  else
      if ( params.pdf_code != 0l )
        LHAPDF::initPDFSet( (int)params.pdf_code, params.pdf_member );
      else
        LHAPDF::initPDFSet( params.pdf_set, LHAPDF::LHGRID, params.pdf_member );
      lhapdf_version = LHAPDF::getVersion();
      pdf_description = LHAPDF::getDescription();
#  endif
      replace_all( pdf_description, ". ", ".\n  " );
      CG_INFO( "Partonic" ) << "Partonic structure functions evaluator successfully built.\n"
        << " * LHAPDF version: " << lhapdf_version << "\n"
        << " * number of flavours: " << params.num_flavours << "\n"
        << " * PDF set: " << params.pdf_set << "\n"
        << ( pdf_description.empty() ? "" : "  "+pdf_description+"\n" )
        << " * PDF member: " << params.pdf_member << ( pdf_type.empty() ? "" : " ("+pdf_type+")" ) << "\n"
        << " * quarks mode: " << params.mode;
      initialised_ = true;
#else
      throw CG_FATAL( "Partonic" ) << "LHAPDF is not liked to this instance!";
#endif
    }

    Partonic&
    Partonic::operator()( double xbj, double q2 )
    {
#ifdef LIBLHAPDF
      std::pair<double,double> nv = { xbj, q2 };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      F2 = 0.;
      if ( params.num_flavours == 0 || params.num_flavours > 6 )
        return *this;

      if ( !initialised_ )
        initialise();
#  if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION >= 6
      auto& member = *pdfs_[params.pdf_member];
      if ( !member.inPhysicalRangeXQ2( xbj, q2 ) ) {
        CG_WARNING( "Partonic" ) << "(x=" << xbj << ", Q²=" << q2 << " GeV²) "
          << "not in physical range for PDF member " << params.pdf_member << ":\n\t"
          << "  min: (x=" << member.xMin() << ", Q²=" << member.q2Min() << "),\n\t"
          << "  max: (x=" << member.xMax() << ", Q²=" << member.q2Max() << ").";
        return *this;
      }
#  else
      if ( q2 < LHAPDF::getQ2min( params.pdf_member ) || q2 > LHAPDF::getQ2max( params.pdf_member )
        || xbj < LHAPDF::getXmin( params.pdf_member ) || xbj > LHAPDF::getXmax( params.pdf_member ) ) {
        CG_WARNING( "Partonic" ) << "(x=" << xbj << "/Q²=" << q2 << " GeV²) "
          << "not in physical range for PDF member " << params.pdf_member << ":\n"
          << "  min: (x=" << LHAPDF::getXmin( params.pdf_member ) << "/Q²=" << LHAPDF::getQ2min( params.pdf_member ) << "),\n"
          << "  max: (x=" << LHAPDF::getXmax( params.pdf_member ) << "/Q²=" << LHAPDF::getQ2max( params.pdf_member ) << ").";
        return *this;
      }
      const double q = sqrt( q2 );
#  endif

      for ( int i = 0; i < params.num_flavours; ++i ) {
        const double prefactor = 1./9.*Q_TIMES_3[i]*Q_TIMES_3[i];
#  if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION >= 6
        if ( !pdfs_[params.pdf_member]->hasFlavor( QUARK_PDGS[i] ) )
          throw CG_FATAL( "Partonic" ) << "Flavour " << QUARK_PDGS[i] << " is unsupported!";
        const double xq = member.xfxQ2( QUARK_PDGS[i], xbj, q2 );
        const double xqbar = member.xfxQ2( -QUARK_PDGS[i], xbj, q2 );
#  else
        const double xq = LHAPDF::xfx( xbj, q, QUARK_PDGS[i] );
        const double xqbar = LHAPDF::xfx( xbj, q, -QUARK_PDGS[i] );
#  endif
        switch ( params.mode ) {
          case Parameters::Mode::full:
            F2 += prefactor*( xq+xqbar ); break;
          case Parameters::Mode::valence:
            F2 += prefactor*( xq-xqbar ); break;
          case Parameters::Mode::sea:
            F2 += prefactor*( 2.*xqbar ); break;
        }
      }
#else
      throw CG_FATAL( "Partonic" ) << "LHAPDF is not liked to this instance!";
#endif

      return *this;
    }
  }

  std::ostream&
  operator<<( std::ostream& os, const strfun::Partonic::Parameters::Mode& mode )
  {
    switch ( mode ) {
      case strfun::Partonic::Parameters::Mode::full: return os << "all quarks";
      case strfun::Partonic::Parameters::Mode::valence: return os << "valence quarks";
      case strfun::Partonic::Parameters::Mode::sea: return os << "sea quarks";
    }
    return os;
  }
}
