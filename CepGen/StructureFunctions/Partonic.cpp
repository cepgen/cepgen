#include "CepGen/StructureFunctions/Partonic.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#ifdef LIBLHAPDF
#  if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
#    define LHAPDF_GE_6 1
#  endif
#endif

namespace cepgen
{
  namespace strfun
  {
    constexpr std::array<short,6> Partonic::Q_TIMES_3, Partonic::QUARK_PDGS;

    Partonic::Partonic( const ParametersList& params ) :
      Parameterisation( params ),
      pdf_set_     ( params.get<std::string>( "pdfSet", "cteq6" ) ),
      num_flavours_( params.get<int>( "numFlavours", 4 ) ),
      pdf_code_    ( params.get<int>( "pdfCode", 0 ) ),
      pdf_member_  ( params.get<int>( "pdfMember", 0 ) ),
      mode_        ( (Mode)params.get<int>( "mode", (int)Mode::full ) ),
      initialised_( false )
    {}

    Partonic::Partonic( const char* set, unsigned short member, const Mode& mode ) :
      Parameterisation( ParametersList().set<int>( "id", (int)Type::Partonic ) ),
      pdf_set_( set ), num_flavours_( 4 ), pdf_code_( 0 ), pdf_member_( member ), mode_( mode ),
      initialised_( false )
    {}

    std::string
    Partonic::description() const
    {
      std::ostringstream os;
      os << "Partonic{" << pdf_set_ << ",m=" << pdf_member_ << ",mode=" << mode_ << "}";
      return os.str();
    }

    void
    Partonic::initialise()
    {
      if ( initialised_ )
        return;
#ifdef LIBLHAPDF
      std::string lhapdf_version, pdf_description, pdf_type;
#  ifdef LHAPDF_GE_6
      try {
        //--- check if PDF code is set
        if ( pdf_code_ != 0l ) {
          auto pdf = LHAPDF::lookupPDF( pdf_code_ );
          if ( pdf.second != 0 )
            throw CG_FATAL( "Partonic" ) << "Failed to retrieve PDFset with id=" << pdf_code_ << "!";
          if ( !pdf_set_.empty() && pdf_set_ != pdf.first )
            CG_WARNING( "Partonic" ) << "PDF set name changed from \"" << pdf_set_ << "\" to \"" << pdf.first << "\".";
          pdf_set_ = pdf.first;
        }
        lha_pdf_set_ = LHAPDF::PDFSet( pdf_set_ );
        lha_pdf_set_.mkPDFs<std::unique_ptr<LHAPDF::PDF> >( pdfs_ );
        lhapdf_version = LHAPDF::version();
        pdf_description = lha_pdf_set_.description();
        pdf_type = pdfs_[pdf_member_]->type();
      } catch ( const LHAPDF::Exception& e ) {
        throw CG_FATAL( "Partonic" )
          << "Caught LHAPDF exception:\n\t"
          << e.what();
      }
#  else
      if ( pdf_code_ != 0l )
        LHAPDF::initPDFSet( (int)pdf_code_, pdf_member_ );
      else
        LHAPDF::initPDFSet( pdf_set_, LHAPDF::LHGRID, pdf_member_ );
      lhapdf_version = LHAPDF::getVersion();
#  endif
      replace_all( pdf_description, ". ", ".\n  " );
      CG_INFO( "Partonic" ) << "Partonic structure functions evaluator successfully built.\n"
        << " * LHAPDF version: " << lhapdf_version << "\n"
        << " * number of flavours: " << num_flavours_ << "\n"
        << " * quarks mode: " << mode_ << "\n"
        << " * PDF set: " << pdf_set_ << "\n"
        << " * PDF member: " << pdf_member_ << ( pdf_type.empty() ? "" : " ("+pdf_type+")" ) << "\n"
#  ifdef LHAPDF_GE_6
        << ( pdf_description.empty() ? "" : "  "+pdf_description );
#  else
      ; LHAPDF::getDescription();
#  endif
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
      if ( num_flavours_ == 0 || num_flavours_ > 6 )
        return *this;

      if ( !initialised_ )
        initialise();
#  ifdef LHAPDF_GE_6
      auto& member = *pdfs_[pdf_member_];
      if ( !member.inPhysicalRangeXQ2( xbj, q2 ) ) {
        CG_WARNING( "Partonic" ) << "(x=" << xbj << ", Q²=" << q2 << " GeV²) "
          << "not in physical range for PDF member " << pdf_member_ << ":\n\t"
          << "  min: (x=" << member.xMin() << ", Q²=" << member.q2Min() << "),\n\t"
          << "  max: (x=" << member.xMax() << ", Q²=" << member.q2Max() << ").";
        return *this;
      }
#  else
      if ( q2 < LHAPDF::getQ2min( pdf_member_ ) || q2 > LHAPDF::getQ2max( pdf_member_ )
        || xbj < LHAPDF::getXmin( pdf_member_ ) || xbj > LHAPDF::getXmax( pdf_member_ ) ) {
        CG_WARNING( "Partonic" ) << "(x=" << xbj << "/Q²=" << q2 << " GeV²) "
          << "not in physical range for PDF member " << pdf_member_ << ":\n"
          << "  min: (x=" << LHAPDF::getXmin( pdf_member_ ) << "/Q²=" << LHAPDF::getQ2min( pdf_member_ ) << "),\n"
          << "  max: (x=" << LHAPDF::getXmax( pdf_member_ ) << "/Q²=" << LHAPDF::getQ2max( pdf_member_ ) << ").";
        return *this;
      }
      const double q = sqrt( q2 );
#  endif

      for ( int i = 0; i < num_flavours_; ++i ) {
        const double prefactor = 1./9.*Q_TIMES_3[i]*Q_TIMES_3[i];
#  ifdef LHAPDF_GE_6
        if ( !pdfs_[pdf_member_]->hasFlavor( QUARK_PDGS[i] ) )
          throw CG_FATAL( "Partonic" ) << "Flavour " << QUARK_PDGS[i] << " is unsupported!";
        const double xq = member.xfxQ2( QUARK_PDGS[i], xbj, q2 );
        const double xqbar = member.xfxQ2( -QUARK_PDGS[i], xbj, q2 );
#  else
        const double xq = LHAPDF::xfx( xbj, q, QUARK_PDGS[i] );
        const double xqbar = LHAPDF::xfx( xbj, q, -QUARK_PDGS[i] );
#  endif
        switch ( mode_ ) {
          case Mode::full:
            F2 += prefactor*( xq+xqbar ); break;
          case Mode::valence:
            F2 += prefactor*( xq-xqbar ); break;
          case Mode::sea:
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
  operator<<( std::ostream& os, const strfun::Partonic::Mode& mode )
  {
    switch ( mode ) {
      case strfun::Partonic::Mode::full: return os << "all quarks";
      case strfun::Partonic::Mode::valence: return os << "valence quarks";
      case strfun::Partonic::Mode::sea: return os << "sea quarks";
    }
    return os;
  }
}

#ifdef LHAPDF_GE_6
#undef LHAPDF_GE_6
#endif

