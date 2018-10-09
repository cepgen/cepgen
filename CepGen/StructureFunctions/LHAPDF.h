#ifndef CepGen_StructureFunctions_LHAPDF_h
#define CepGen_StructureFunctions_LHAPDF_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

#ifdef LIBLHAPDF
#include "LHAPDF/LHAPDF.h"
#endif

#include <array>

namespace CepGen
{
  namespace SF
  {
    /// Generic partonic level perturbative structure functions built from an external PDFs grid
    class LHAPDF : public Parameterisation
    {
      public:
        /// Model parameterisation
        struct Parameters
        {
          /// Standard (usual CTEQ6) constructor
          Parameters();
          /// Number of quark flavours considered in the SF building
          unsigned short num_flavours;
          /// String-type PDF identifier (default)
          std::string pdf_set;
          /// Integer-type PDF identifier (if no string version is provided)
          unsigned long pdf_code;
          /// PDF set used
          unsigned short pdf_member;
          /// Quarks types
          enum class Mode { full = 0, valence = 1, sea = 2 };
          /// Quarks types considered in the SF building
          Mode mode;
        };
        /// Build a calculator from its Parameters object
        explicit LHAPDF( const Parameters& param = Parameters() );
        /// Build a calculator from a set, its member, and the contributing quarks
        explicit LHAPDF( const char* set, unsigned short member = 0, const Parameters::Mode& mode = Parameters::Mode::full );
        LHAPDF& operator()( double xbj, double q2 ) override;
        /// Parameterisation used in this SFs calculator
        Parameters params;

      private:
        std::string description() const override;
        void initialise();
        bool initialised_;

#ifdef LIBLHAPDF
#  if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION >= 6
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
  std::ostream& operator<<( std::ostream& os, const SF::LHAPDF::Parameters::Mode& mode );
}

#endif
