#ifndef CepGen_StructureFunctions_Partonic_h
#define CepGen_StructureFunctions_Partonic_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

#ifdef LIBLHAPDF
#include "LHAPDF/LHAPDF.h"
#endif

#include <array>

namespace cepgen
{
  namespace strfun
  {
    /// Generic partonic level perturbative structure functions built from an external PDFs grid
    class Partonic : public Parameterisation
    {
      public:
        /// Quarks types
        enum class Mode { full = 0, valence = 1, sea = 2 };
        /// Build a calculator from its Parameters object
        explicit Partonic( const ParametersList& params = ParametersList() );
        /// Build a calculator from a set, its member, and the contributing quarks
        explicit Partonic( const char* set, unsigned short member = 0, const Mode& mode = Mode::full );
        Partonic& operator()( double xbj, double q2 ) override;

      private:
        std::string description() const override;
        void initialise();
        /// String-type PDF identifier (default)
        std::string pdf_set_;
        /// Number of quark flavours considered in the SF building
        unsigned short num_flavours_;
        /// Integer-type PDF identifier (if no string version is provided)
        unsigned long pdf_code_;
        /// PDF set used
        unsigned short pdf_member_;
        /// Quarks types considered in the SF building
        Mode mode_;
        bool initialised_;

#ifdef LIBLHAPDF
#  if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION >= 6
        LHAPDF::PDFSet lha_pdf_set_;
        std::vector<std::unique_ptr<LHAPDF::PDF> > pdfs_;
#  endif
#endif
        static constexpr std::array<short,6> QUARK_PDGS = { { 1, 2, 3, 4, 5, 6 } };
        static constexpr std::array<short,6> Q_TIMES_3 = { {
          -1 /*d*/, 2 /*u*/,
          -1 /*s*/, 2 /*c*/,
          -1 /*b*/, 2 /*t*/
        } };
    };
  }
  std::ostream& operator<<( std::ostream& os, const strfun::Partonic::Mode& mode );
}

#endif

