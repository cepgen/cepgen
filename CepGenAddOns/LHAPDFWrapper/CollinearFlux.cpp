/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen::lhapdf {
  class CollinearFlux final : public cepgen::CollinearFlux {
  public:
    explicit CollinearFlux(const ParametersList& params)
        : cepgen::CollinearFlux(params),
          pdf_(LHAPDF::mkPDF(steer<std::string>("set"), steer<int>("member"))),
          pdgid_(steerAs<int, pdgid_t>("partonPdgId")),
          extrapolate_pdf_(steer<bool>("extrapolatePDF")) {
      const auto& pdf_set = steer<std::string>("set");
      if (!pdf_)
        throw CG_FATAL("lhapdf:CollinearFlux") << "Failed to initialise the LHAPDF evaluator!\n"
                                               << "Parameters: " << params_;
      if (extrapolate_pdf_ && pdf_->hasFlavor(pdgid_))
        CG_WARNING("lhapdf:CollinearFlux") << "Asked to retrieve distribution from sum imbalance of other "
                                              "contributions although the distribution is present in the '"
                                           << pdf_set << "' PDF set.\n\t"
                                           << "You may want to steer the 'extrapolatePDF' parameter to 'false'?";
      if (!extrapolate_pdf_ && !pdf_->hasFlavor(pdgid_))
        throw CG_FATAL("lhapdf:CollinearFlux")
            << "PDF set '" << pdf_set << "' does not contain parton with PDG identifier=" << pdgid_ << "!\n"
            << "PDGs handled: " << pdf_->flavors() << ".";

      CG_INFO("lhapdf:CollinearFlux") << "LHAPDF evaluator for collinear parton (" << (PDG::Id)pdgid_
                                      << ") flux initialised.\n\t"
                                      << "PDF set: " << steer<std::string>("set") << " (flavours: " << pdf_->flavors()
                                      << "), member: " << steer<int>("member") << ".\n\t"
                                      << "x range: " << Limits{pdf_->xMin(), pdf_->xMax()} << ", "
                                      << "Q^2 range: " << Limits{pdf_->q2Min(), pdf_->q2Max()} << " GeV^2.\n\t"
                                      << "Extrapolated from other flavours? " << extrapolate_pdf_ << ".";
    }

    static ParametersDescription description() {
      auto desc = cepgen::CollinearFlux::description();
      desc.setDescription("LHAPDF coll.flux");
      desc.add<std::string>("set", "LUXqed17_plus_PDF4LHC15_nnlo_100").setDescription("PDFset to use");
      desc.add<int>("member", 0).setDescription("PDF member");
      desc.addAs<int, pdgid_t>("partonPdgId", PDG::photon).setDescription("parton PDG identifier");
      desc.add<bool>("extrapolatePDF", false)
          .setDescription("has the PDF? or extrapolate distribution from sum imbalance of other contributions?");
      return desc;
    }

    pdgid_t partonPdgId() const override { return pdgid_; }
    bool fragmenting() const override { return true; }
    double mass2() const override { return mp2_; }

    double fluxQ2(double x, double q2) const override {
      if (x == 0. || !pdf_->inPhysicalRangeXQ2(x, q2))
        return 0.;
      if (!extrapolate_pdf_)  // has parton PDF
        return pdf_->xfxQ2((int)pdgid_, x, q2);
      // extrapolate from other flavours imbalance
      double xf = 1.;
      for (const auto& flav : pdf_->xfxQ2(x, q2))
        if (flav.first != (int)pdgid_)
          xf -= flav.second;
      return xf;
    }

  private:
    const std::unique_ptr<LHAPDF::PDF> pdf_;
    const pdgid_t pdgid_;
    const bool extrapolate_pdf_;
  };
}  // namespace cepgen::lhapdf
using LHAPDFCollinearFlux = cepgen::lhapdf::CollinearFlux;
REGISTER_COLLINEAR_FLUX("lhapdf", LHAPDFCollinearFlux);
