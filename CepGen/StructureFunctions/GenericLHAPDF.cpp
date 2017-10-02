#include "GenericLHAPDF.h"

namespace CepGen
{
  namespace SF
  {
    GenericLHAPDF::GenericLHAPDF( const char* set )
    {
      initialise( set );
    }

    void
    GenericLHAPDF::initialise( const char* set )
    {
#if LHAPDF_MAJOR_VERSION==6
      pdf_set_ = LHAPDF::PDFSet( set );
      pdfs_ = pdf_set_.mkPDFs();
#else
      LHAPDF::initPDFSet( set, LHAPDF::LHGRID, 0 );
#endif
    }

    StructureFunctions
    GenericLHAPDF::operator()( double q2, double xbj ) const
    {
      StructureFunctions pdf;

      //const LHAPDF::PDFSet set( "MRST2004qed_proton" );
      std::array<double,6> qtimes3 = { -1.0 /*d*/, 2.0 /*u*/, -1.0 /*s*/, 2.0 /*c*/, -1.0 /*b*/, 2.0 /*t*/ };

      double sf = 0.;
      //if ( q2 < 1.69 ) return pdf;

      for ( int i = 0; i < 4; ++i ) {
        double xq = 0., xqbar = 0.;
#if LHAPDF_MAJOR_VERSION==6
        xq = pdfs_[0]->xfxQ2( i, xbj, q2 );
        xqbar = pdfs_[0]->xfxQ2( -i, xbj, q2 );
#else
        xq = LHAPDF::xfx( xbj, q2, i+1 );
        xqbar = LHAPDF::xfx( xbj, q2, -i-1 );
#endif
        sf += qtimes3[i]*qtimes3[i]/9. * ( xq + xqbar );
      }

      pdf.F2 = sf;
      return pdf;
    }
  }
}
