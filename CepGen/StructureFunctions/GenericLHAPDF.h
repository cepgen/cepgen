#ifndef CepGen_StructureFunctions_CTEQ_h
#define CepGen_StructureFunctions_CTEQ_h

#include "StructureFunctions.h"

#ifdef LIBLHAPDF
#include "LHAPDF/LHAPDF.h"
#endif

#include <array>

namespace CepGen
{
  namespace SF
  {
    /// Generic, tree-level import of structure functions from an external PDFs grid
    class GenericLHAPDF
    {
      public:
        GenericLHAPDF( const char* set );
        StructureFunctions operator()( double q2, double xbj, unsigned short num_flavours = 4 ) const;

      private:
        void initialise( const char* set );

#ifdef LIBLHAPDF
#if LHAPDF_MAJOR_VERSION==6
        LHAPDF::PDFSet pdf_set_;
        std::vector<LHAPDF::PDF*> pdfs_;
#endif
#endif
        static constexpr std::array<double,6> qtimes3_ = { {
          -1.0 /*d*/, 2.0 /*u*/,
          -1.0 /*s*/, 2.0 /*c*/,
          -1.0 /*b*/, 2.0 /*t*/
        } };
    };
  }
}

#endif
