#include "GenericLHAPDF.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  namespace SF
  {
    constexpr std::array<double,6> GenericLHAPDF::qtimes3_;

    GenericLHAPDF::GenericLHAPDF() :
      StructureFunctions( Type::GenericLHAPDF )
    {}

    GenericLHAPDF::GenericLHAPDF( const char* set ) :
      StructureFunctions( Type::GenericLHAPDF )
    {
      initialise( set );
    }

    void
    GenericLHAPDF::initialise( const char* set )
    {
#ifdef LIBLHAPDF
#  if LHAPDF_MAJOR_VERSION == 6
      pdf_set_ = LHAPDF::PDFSet( set );
      pdfs_ = pdf_set_.mkPDFs();
#  else
      LHAPDF::initPDFSet( set, LHAPDF::LHGRID, 0 );
#  endif
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
      if ( num_flavours == 0 || num_flavours > 6 )
        return *this;

#ifdef LIBLHAPDF
      for ( int i = 0; i < num_flavours; ++i ) {
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
