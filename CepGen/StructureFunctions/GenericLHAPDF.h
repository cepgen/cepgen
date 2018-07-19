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
    class GenericLHAPDF : public StructureFunctions
    {
      public:
        struct Parameterisation
        {
          Parameterisation();
          static Parameterisation cteq6();
          unsigned short num_flavours;
          std::string pdf_set;
          unsigned short pdf_member;
        };

        explicit GenericLHAPDF( const Parameterisation& param = Parameterisation::cteq6() );
        explicit GenericLHAPDF( const char* set );
        GenericLHAPDF& operator()( double q2, double xbj ) override;

        Parameterisation params;

      private:
        void initialise();
        bool initialised_;

#ifdef LIBLHAPDF
#  if LHAPDF_MAJOR_VERSION == 6
        LHAPDF::PDFSet pdf_set_;
        std::vector<std::unique_ptr<LHAPDF::PDF> > pdfs_;
#  endif
#endif
        static constexpr std::array<short,6> qtimes3_ = { {
          -1 /*d*/, 2 /*u*/,
          -1 /*s*/, 2 /*c*/,
          -1 /*b*/, 2 /*t*/
        } };
    };
  }
}

#endif
