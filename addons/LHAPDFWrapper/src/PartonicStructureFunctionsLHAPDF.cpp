/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/PartonicParameterisation.h"
#include "CepGen/Utils/String.h"

#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
#define LHAPDF_GE_6 1
#endif

using namespace std::string_literals;

namespace cepgen::strfun {
  /// Generic partonic level perturbative structure functions built from an external PDFs grid
  class LHAPDFPartonic : public PartonicParameterisation {
  public:
    /// Quarks types
    enum class Mode { full = 0, valence = 1, sea = 2 };
    /// Build a calculator from its Parameters object
    explicit LHAPDFPartonic(const ParametersList& params)
        : PartonicParameterisation(params),
          pdf_set_(steer<std::string>("pdfSet")),
          pdf_code_(steer<int>("pdfCode")),
          pdf_member_(steer<int>("pdfMember")) {}

    static ParametersDescription description() {
      auto desc = PartonicParameterisation::description();
      desc.setDescription("LHAPDF (partonic)");
      desc.add("pdfSet", "cteq66"s).setDescription("PDF modelling to be considered");
      desc.add("pdfCode", 0);
      desc.add("pdfMember", 0);
      return desc;
    }

  private:
    void initialise();
    double evalxQ2(int flavour, double xbj, double q2) override;
    /// String-type PDF identifier (default)
    std::string pdf_set_;
    /// Integer-type PDF identifier (if no string version is provided)
    int pdf_code_;
    /// PDF set used
    int pdf_member_;
    bool initialised_{false};

#ifdef LHAPDF_GE_6
    LHAPDF::PDFSet lha_pdf_set_;
    std::vector<std::unique_ptr<LHAPDF::PDF> > pdfs_;
#endif
  };

  void LHAPDFPartonic::initialise() {
    if (initialised_)
      return;
    std::string lhapdf_version, pdf_description, pdf_type;
#ifdef LHAPDF_GE_6
    try {
      //--- check if PDF code is set
      if (pdf_code_ != 0) {
        auto pdf = LHAPDF::lookupPDF(pdf_code_);
        if (pdf.second != 0)
          throw CG_FATAL("LHAPDFPartonic") << "Failed to retrieve PDFset with id=" << pdf_code_ << "!";
        if (!pdf_set_.empty() && pdf_set_ != pdf.first)
          CG_WARNING("LHAPDFPartonic") << "PDF set name changed from \"" << pdf_set_ << "\" to \"" << pdf.first
                                       << "\".";
        pdf_set_ = pdf.first;
      }
      lha_pdf_set_ = LHAPDF::PDFSet(pdf_set_);
      lha_pdf_set_.mkPDFs<std::unique_ptr<LHAPDF::PDF> >(pdfs_);
      lhapdf_version = LHAPDF::version();
      pdf_description = utils::replaceAll(lha_pdf_set_.description(), {{"\\n", "\n"}, {". ", ".\n  "}});
      pdf_type = pdfs_[pdf_member_]->type();
    } catch (const LHAPDF::Exception& e) {
      throw CG_FATAL("LHAPDFPartonic") << "Caught LHAPDF exception:\n\t" << e.what();
    }
#else
    if (pdf_code_ != 0)
      LHAPDF::initPDFSet(pdf_code_, pdf_member_);
    else
      LHAPDF::initPDFSet(pdf_set_, LHAPDF::LHGRID, pdf_member_);
    lhapdf_version = LHAPDF::getVersion();
#endif
    CG_INFO("LHAPDFPartonic") << "Partonic structure functions evaluator successfully built.\n"
                              << " * LHAPDF version: " << lhapdf_version << "\n"
                              << " * number of flavours: " << num_flavours_ << "\n"
                              << " * quarks mode: " << mode_ << "\n"
                              << " * PDF set: " << pdf_set_ << "\n"
                              << " * PDF member: " << pdf_member_ << (pdf_type.empty() ? "" : " (" + pdf_type + ")")
                              << "\n"
                              << (pdf_description.empty() ? "" : "  " + pdf_description);
#ifndef LHAPDF_GE_6
    LHAPDF::getDescription();
#endif
    initialised_ = true;
  }

  double LHAPDFPartonic::evalxQ2(int flavour, double xbj, double q2) {
    if (!initialised_)
      initialise();
#ifdef LHAPDF_GE_6
    auto& member = *pdfs_[pdf_member_];
    if (!member.inPhysicalRangeXQ2(xbj, q2)) {
      CG_WARNING("LHAPDFPartonic") << "(x=" << xbj << ", Q²=" << q2 << " GeV²) "
                                   << "not in physical range for PDF member " << pdf_member_ << ":\n\t"
                                   << "  min: (x=" << member.xMin() << ", Q²=" << member.q2Min() << "),\n\t"
                                   << "  max: (x=" << member.xMax() << ", Q²=" << member.q2Max() << ").";
      return 0.;
    }
    if (!pdfs_[pdf_member_]->hasFlavor(flavour))
      throw CG_FATAL("LHAPDFPartonic") << "Flavour " << flavour << " is unsupported!";
    return member.xfxQ2(flavour, xbj, q2);
#else
    if (q2 < LHAPDF::getQ2min(pdf_member_) || q2 > LHAPDF::getQ2max(pdf_member_) ||
        xbj < LHAPDF::getXmin(pdf_member_) || xbj > LHAPDF::getXmax(pdf_member_)) {
      CG_WARNING("LHAPDFPartonic") << "(x=" << xbj << "/Q²=" << q2 << " GeV²) "
                                   << "not in physical range for PDF member " << pdf_member_ << ":\n"
                                   << "  min: (x=" << LHAPDF::getXmin(pdf_member_)
                                   << "/Q²=" << LHAPDF::getQ2min(pdf_member_) << "),\n"
                                   << "  max: (x=" << LHAPDF::getXmax(pdf_member_)
                                   << "/Q²=" << LHAPDF::getQ2max(pdf_member_) << ").";
      return 0.;
    }
    return LHAPDF::xfx(xbj, std::sqrt(q2), flavour);
#endif
  }
}  // namespace cepgen::strfun

#ifdef LHAPDF_GE_6
#undef LHAPDF_GE_6
#endif

using cepgen::strfun::LHAPDFPartonic;
REGISTER_STRFUN("lhapdf", 401, LHAPDFPartonic);
