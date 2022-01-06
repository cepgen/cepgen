/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <LHAPDF/LHAPDF.h>

#include <array>
#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/String.h"

#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
#define LHAPDF_GE_6 1
#endif

namespace cepgen {
  namespace strfun {
    /// Generic partonic level perturbative structure functions built from an external PDFs grid
    class Partonic : public Parameterisation {
    public:
      /// Quarks types
      enum class Mode { full = 0, valence = 1, sea = 2 };
      /// Build a calculator from its Parameters object
      explicit Partonic(const ParametersList& params = ParametersList());
      /// Build a calculator from a set, its member, and the contributing quarks
      explicit Partonic(const char* set, unsigned short member = 0, const Mode& mode = Mode::full);

      static ParametersDescription description();

      Partonic& eval(double xbj, double q2) override;
      std::string describe() const override;

    private:
      void initialise();
      /// String-type PDF identifier (default)
      std::string pdf_set_;
      /// Number of quark flavours considered in the SF building
      unsigned short num_flavours_;
      /// Integer-type PDF identifier (if no string version is provided)
      int pdf_code_;
      /// PDF set used
      int pdf_member_;
      /// Quarks types considered in the SF building
      Mode mode_;
      bool initialised_{false};

#ifdef LHAPDF_GE_6
      LHAPDF::PDFSet lha_pdf_set_;
      std::vector<std::unique_ptr<LHAPDF::PDF> > pdfs_;
#endif
      static constexpr std::array<short, 6> QUARK_PDGS = {{1, 2, 3, 4, 5, 6}};
      static constexpr std::array<short, 6> Q_TIMES_3 = {{
          -1 /*d*/, 2 /*u*/, -1 /*s*/, 2 /*c*/, -1 /*b*/, 2 /*t*/
      }};
    };

    constexpr std::array<short, 6> Partonic::Q_TIMES_3, Partonic::QUARK_PDGS;

    std::ostream& operator<<(std::ostream& os, const strfun::Partonic::Mode& mode) {
      switch (mode) {
        case strfun::Partonic::Mode::full:
          return os << "all quarks";
        case strfun::Partonic::Mode::valence:
          return os << "valence quarks";
        case strfun::Partonic::Mode::sea:
          return os << "sea quarks";
      }
      return os;
    }

    Partonic::Partonic(const ParametersList& params)
        : Parameterisation(params),
          pdf_set_(params.get<std::string>("pdfSet")),
          num_flavours_(params.get<int>("numFlavours")),
          pdf_code_(params.get<int>("pdfCode")),
          pdf_member_(params.get<int>("pdfMember")),
          mode_(params.getAs<int, Mode>("mode")) {}

    Partonic::Partonic(const char* set, unsigned short member, const Mode& mode)
        : Parameterisation(ParametersList().setName<int>((int)Type::Partonic)),
          pdf_set_(set),
          num_flavours_(4),
          pdf_code_(0),
          pdf_member_(member),
          mode_(mode),
          initialised_(false) {}

    std::string Partonic::describe() const {
      std::ostringstream os;
      os << "Partonic{" << pdf_set_ << ",m=" << pdf_member_ << ",mode=" << mode_ << "}";
      return os.str();
    }

    void Partonic::initialise() {
      if (initialised_)
        return;
      std::string lhapdf_version, pdf_description, pdf_type;
#ifdef LHAPDF_GE_6
      try {
        //--- check if PDF code is set
        if (pdf_code_ != 0) {
          auto pdf = LHAPDF::lookupPDF(pdf_code_);
          if (pdf.second != 0)
            throw CG_FATAL("Partonic") << "Failed to retrieve PDFset with id=" << pdf_code_ << "!";
          if (!pdf_set_.empty() && pdf_set_ != pdf.first)
            CG_WARNING("Partonic") << "PDF set name changed from \"" << pdf_set_ << "\" to \"" << pdf.first << "\".";
          pdf_set_ = pdf.first;
        }
        lha_pdf_set_ = LHAPDF::PDFSet(pdf_set_);
        lha_pdf_set_.mkPDFs<std::unique_ptr<LHAPDF::PDF> >(pdfs_);
        lhapdf_version = LHAPDF::version();
        pdf_description = utils::replace_all(lha_pdf_set_.description(), ". ", ".\n  ");
        ;
        pdf_type = pdfs_[pdf_member_]->type();
      } catch (const LHAPDF::Exception& e) {
        throw CG_FATAL("Partonic") << "Caught LHAPDF exception:\n\t" << e.what();
      }
#else
      if (pdf_code_ != 0)
        LHAPDF::initPDFSet(pdf_code_, pdf_member_);
      else
        LHAPDF::initPDFSet(pdf_set_, LHAPDF::LHGRID, pdf_member_);
      lhapdf_version = LHAPDF::getVersion();
#endif
      CG_INFO("Partonic") << "Partonic structure functions evaluator successfully built.\n"
                          << " * LHAPDF version: " << lhapdf_version << "\n"
                          << " * number of flavours: " << num_flavours_ << "\n"
                          << " * quarks mode: " << mode_ << "\n"
                          << " * PDF set: " << pdf_set_ << "\n"
                          << " * PDF member: " << pdf_member_ << (pdf_type.empty() ? "" : " (" + pdf_type + ")") << "\n"
                          << (pdf_description.empty() ? "" : "  " + pdf_description);
#ifndef LHAPDF_GE_6
      LHAPDF::getDescription();
#endif
      initialised_ = true;
    }

    Partonic& Partonic::eval(double xbj, double q2) {
      F2 = 0.;
      if (num_flavours_ == 0 || num_flavours_ > QUARK_PDGS.size()) {
        CG_WARNING("Partonic") << "Invalid number of flavours (" << num_flavours_ << " selected.";
        return *this;
      }

      if (!initialised_)
        initialise();
#ifdef LHAPDF_GE_6
      auto& member = *pdfs_[pdf_member_];
      if (!member.inPhysicalRangeXQ2(xbj, q2)) {
        CG_WARNING("Partonic") << "(x=" << xbj << ", Q²=" << q2 << " GeV²) "
                               << "not in physical range for PDF member " << pdf_member_ << ":\n\t"
                               << "  min: (x=" << member.xMin() << ", Q²=" << member.q2Min() << "),\n\t"
                               << "  max: (x=" << member.xMax() << ", Q²=" << member.q2Max() << ").";
        return *this;
      }
#else
      if (q2 < LHAPDF::getQ2min(pdf_member_) || q2 > LHAPDF::getQ2max(pdf_member_) ||
          xbj < LHAPDF::getXmin(pdf_member_) || xbj > LHAPDF::getXmax(pdf_member_)) {
        CG_WARNING("Partonic") << "(x=" << xbj << "/Q²=" << q2 << " GeV²) "
                               << "not in physical range for PDF member " << pdf_member_ << ":\n"
                               << "  min: (x=" << LHAPDF::getXmin(pdf_member_)
                               << "/Q²=" << LHAPDF::getQ2min(pdf_member_) << "),\n"
                               << "  max: (x=" << LHAPDF::getXmax(pdf_member_)
                               << "/Q²=" << LHAPDF::getQ2max(pdf_member_) << ").";
        return *this;
      }
      const double q = sqrt(q2);
#endif

      for (int i = 0; i < num_flavours_; ++i) {
        const double prefactor = 1. / 9. * Q_TIMES_3.at(i) * Q_TIMES_3.at(i);
#ifdef LHAPDF_GE_6
        if (!pdfs_[pdf_member_]->hasFlavor(QUARK_PDGS.at(i)))
          throw CG_FATAL("Partonic") << "Flavour " << QUARK_PDGS.at(i) << " is unsupported!";
        const double xq = member.xfxQ2(QUARK_PDGS.at(i), xbj, q2);
        const double xqbar = member.xfxQ2(-QUARK_PDGS.at(i), xbj, q2);
#else
        const double xq = LHAPDF::xfx(xbj, q, QUARK_PDGS.at(i));
        const double xqbar = LHAPDF::xfx(xbj, q, -QUARK_PDGS.at(i));
#endif
        switch (mode_) {
          case Mode::full:
            F2 += prefactor * (xq + xqbar);
            break;
          case Mode::valence:
            F2 += prefactor * (xq - xqbar);
            break;
          case Mode::sea:
            F2 += prefactor * (2. * xqbar);
            break;
        }
      }
      return *this;
    }

    ParametersDescription Partonic::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Partonic structure functions");
      desc.add<std::string>("pdfSet", "cteq6").setDescription("PDF modelling to be considered");
      desc.add<int>("numFlavours", 4).setDescription("Number of parton flavours to consider in summation");
      desc.add<int>("pdfCode", 0);
      desc.add<int>("pdfMember", 0);
      desc.add<int>("mode", (int)Mode::full);
      return desc;
    }
  }  // namespace strfun
}  // namespace cepgen

#ifdef LHAPDF_GE_6
#undef LHAPDF_GE_6
#endif

REGISTER_STRFUN(strfun::Type::Partonic, Partonic, strfun::Partonic)
