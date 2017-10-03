#ifndef CepGen_StructureFunctions_CTEQ_h
#define CepGen_StructureFunctions_CTEQ_h

#include "StructureFunctions.h"

#include "LHAPDF/LHAPDF.h"
#include <array>

namespace CepGen
{
  namespace SF
  {
    class GenericLHAPDF
    {
      public:
        GenericLHAPDF( const char* set );
        StructureFunctions operator()( double q2, double xbj, unsigned short num_flavours = 4 ) const;

      private:
        void initialise( const char* set );

#if LHAPDF_MAJOR_VERSION==6
        LHAPDF::PDFSet pdf_set_;
        std::vector<LHAPDF::PDF*> pdfs_;
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
