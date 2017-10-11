#include "GenericLHAPDF.h"

namespace CepGen
{
  namespace SF
  {
    constexpr std::array<double,6> GenericLHAPDF::qtimes3_;
    GenericLHAPDF::GenericLHAPDF( const char* set )
    {
      initialise( set );
    }

    void
    GenericLHAPDF::initialise( const char* set )
    {
#ifdef LIBLHAPDF
#if LHAPDF_MAJOR_VERSION==6
      pdf_set_ = LHAPDF::PDFSet( set );
      pdfs_ = pdf_set_.mkPDFs();
#else
      LHAPDF::initPDFSet( set, LHAPDF::LHGRID, 0 );
#endif
#else
      FatalError( "LHAPDF is not liked to this instance!" );
#endif
    }

    StructureFunctions
    GenericLHAPDF::operator()( double q2, double xbj, unsigned short num_flavours ) const
    {
      StructureFunctions pdf;

      if ( num_flavours == 0 || num_flavours > 6 ) return pdf;

      //if ( q2 < 1.69 ) return pdf;
#ifdef LIBLHAPDF
      for ( int i = 0; i < num_flavours; ++i ) {
        double xq = 0., xqbar = 0.;
#if LHAPDF_MAJOR_VERSION==6
        xq = pdfs_[0]->xfxQ2( i, xbj, q2 );
        xqbar = pdfs_[0]->xfxQ2( -i, xbj, q2 );
#else
        xq = LHAPDF::xfx( xbj, q2, i+1 );
        xqbar = LHAPDF::xfx( xbj, q2, -i-1 );
#endif
        pdf.F2 += qtimes3_[i]*qtimes3_[i]/9. * ( xq + xqbar );
      }
#else
      FatalError( "LHAPDF is not liked to this instance!" );
#endif

      return pdf;
    }
  }
}
