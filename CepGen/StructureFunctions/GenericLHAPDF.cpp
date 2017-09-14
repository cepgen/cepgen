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
      pdf_set_ = LHAPDF::PDFSet( set );
      pdfs_ = pdf_set_.mkPDFs();      
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
        sf += qtimes3[i]*qtimes3[i]/9. * ( pdfs_[0]->xfxQ2( i, xbj, q2 ) + pdfs_[0]->xfxQ2( -i, xbj, q2 ) );
      }

      pdf.F2 = sf;
      return pdf;
    }
  }
}
