#ifndef CepGen_StructureFunctions_CTEQ_h
#define CepGen_StructureFunctions_CTEQ_h

#include "StructureFunctions.h"
#include "LHAPDF/LHAPDF.h"

namespace CepGen
{
  namespace SF
  {
    class GenericLHAPDF
    {
      public:
        GenericLHAPDF( const char* set );
        StructureFunctions operator()( double q2, double xbj ) const;

      private:
        void initialise( const char* set );

        LHAPDF::PDFSet pdf_set_;
        std::vector<LHAPDF::PDF*> pdfs_;
    };
  }
}

#endif
