#ifndef CepGen_StructureFunctions_LHAPDF_h
#define CepGen_StructureFunctions_LHAPDF_h

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
    class LHAPDF : public StructureFunctions
    {
      public:
        struct Parameterisation
        {
          Parameterisation();
          static Parameterisation cteq6();
          unsigned short num_flavours;
          std::string pdf_set;
          unsigned short pdf_member;
          enum class Mode { full = 0, valence = 1, sea = 2 };
          Mode mode;
        };

        explicit LHAPDF( const Parameterisation& param = Parameterisation::cteq6() );
        explicit LHAPDF( const char* set, unsigned short member = 0, const Parameterisation::Mode& mode = Parameterisation::Mode::full );
        LHAPDF& operator()( double q2, double xbj ) override;

        Parameterisation params;

      private:
        void initialise();
        bool initialised_;

#ifdef LIBLHAPDF
#  if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
        ::LHAPDF::PDFSet pdf_set_;
        std::vector<std::unique_ptr<::LHAPDF::PDF> > pdfs_;
#  endif
#endif
        static constexpr std::array<short,6> pdgid_ = { { 1, 2, 3, 4, 5, 6 } };
        static constexpr std::array<short,6> qtimes3_ = { {
          -1 /*d*/, 2 /*u*/,
          -1 /*s*/, 2 /*c*/,
          -1 /*b*/, 2 /*t*/
        } };
    };
  }
  std::ostream& operator<<( std::ostream& os, const SF::LHAPDF::Parameterisation::Mode& mode );
}

#endif
