/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include <numeric>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/Physics/PDG.h"

using namespace std::string_literals;

namespace cepgen::lhapdf {
  class CollinearFlux final : public cepgen::CollinearFlux {
  public:
    explicit CollinearFlux(const ParametersList& params)
        : cepgen::CollinearFlux(params),
          pdf_(LHAPDF::mkPDF(steer<std::string>("set"), steer<int>("member"))),
          parton_pdgid_(steer<int>("partonPdgId")),
          extrapolate_pdf_(steer<bool>("extrapolatePDF")) {
      const auto& pdf_set = steer<std::string>("set");
      if (!pdf_)
        throw CG_FATAL("lhapdf:CollinearFlux") << "Failed to initialise the LHAPDF evaluator!\n"
                                               << "Parameters: " << params_;
      if (extrapolate_pdf_ && pdf_->hasFlavor(parton_pdgid_))
        CG_WARNING("lhapdf:CollinearFlux") << "Asked to retrieve distribution from sum imbalance of other "
                                              "contributions although the distribution is present in the '"
                                           << pdf_set << "' PDF set.\n\t"
                                           << "You may want to steer the 'extrapolatePDF' parameter to 'false'?";
      if (!extrapolate_pdf_ && !pdf_->hasFlavor(parton_pdgid_))
        throw CG_FATAL("lhapdf:CollinearFlux")
            << "PDF set '" << pdf_set << "' does not contain parton with PDG identifier=" << parton_pdgid_ << "!\n"
            << "PDGs handled: " << pdf_->flavors() << ".";

      CG_INFO("lhapdf:CollinearFlux") << "LHAPDF evaluator for collinear parton ("
                                      << static_cast<PDG::Id>(parton_pdgid_) << ") flux initialised.\n\t"
                                      << "PDF set: " << pdf_set << " (flavours: " << pdf_->flavors()
                                      << "), member: " << pdf_->memberID() << ".\n\t"
                                      << "x range: " << Limits{pdf_->xMin(), pdf_->xMax()} << ", "
                                      << "Q^2 range: " << Limits{pdf_->q2Min(), pdf_->q2Max()} << " GeV^2.\n\t"
                                      << "Extrapolated from other flavours? " << extrapolate_pdf_ << ".";
    }

    static ParametersDescription description() {
      auto desc = cepgen::CollinearFlux::description();
      desc.setDescription("LHAPDF coll.flux");
      desc.add("set", "LUXlep-NNPDF31_nlo_as_0118_luxqed"s).setDescription("PDFset to use");
      desc.add("member", 0).setDescription("PDF member");
      desc.addAs<int>("partonPdgId", PDG::photon).setDescription("parton PDG identifier");
      desc.add("extrapolatePDF", false)
          .setDescription("has the PDF? or extrapolate distribution from sum imbalance of other contributions?");
      return desc;
    }

    spdgid_t partonPdgId() const override { return parton_pdgid_; }
    bool fragmenting() const override { return true; }
    double mass2() const override { return mp2_; }

    double fluxQ2(double x, double q2) const override {
      if (x == 0. || !pdf_->inPhysicalRangeXQ2(x, q2))
        return 0.;
      if (!extrapolate_pdf_)  // has parton PDF
        return pdf_->xfxQ2(parton_pdgid_, x, q2);
      // extrapolate from other flavours imbalance
      double xf = 1.;
      for (const auto& flav : pdf_->xfxQ2(x, q2))
        if (flav.first != parton_pdgid_)
          xf -= flav.second;
      return xf;
    }

  private:
    const std::unique_ptr<LHAPDF::PDF> pdf_;
    const spdgid_t parton_pdgid_;
    const bool extrapolate_pdf_;
  };
}  // namespace cepgen::lhapdf
using LHAPDFCollinearFlux = cepgen::lhapdf::CollinearFlux;
REGISTER_COLLINEAR_FLUX("lhapdf", LHAPDFCollinearFlux);
