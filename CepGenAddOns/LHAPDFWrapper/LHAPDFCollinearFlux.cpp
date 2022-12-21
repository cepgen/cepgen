/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"

namespace cepgen {
  namespace collflux {
    class LHAPDFCollinearFlux : public Parameterisation {
    public:
      explicit LHAPDFCollinearFlux(const ParametersList& params)
          : Parameterisation(params),
            pdf_(LHAPDF::mkPDF(steer<std::string>("set"), steer<int>("member"))),
            pdgid_(steerAs<int, pdgid_t>("pdgId")),
            from_remnant_(steer<bool>("fromRemnant")) {
        const auto& pdf_set = steer<std::string>("set");
        if (!pdf_)
          throw CG_FATAL("LHAPDFCollinearFlux") << "Failed to initialise the LHAPDF evaluator!\n"
                                                << "Parameters: " << params_;
        if (from_remnant_ && pdf_->hasFlavor(pdgid_))
          CG_WARNING("LHAPDFCollinearFlux") << "Asked to retrieve distribution from sum imbalance of other "
                                               "contributions although the distribution is present in the '"
                                            << pdf_set << "' PDF set.";
        if (!from_remnant_ && !pdf_->hasFlavor(pdgid_))
          throw CG_FATAL("LHAPDFCollinearFlux")
              << "PDF set '" << pdf_set << "' does not contain parton with PDG identifier=" << pdgid_ << "!\n"
              << "PDGs handled: " << pdf_->flavors() << ".";

        CG_INFO("LHAPDFCollinearFlux") << "LHAPDF evaluator for collinear parton flux initialised.\n\t"
                                       << "Parton PDG identifier: " << pdgid_ << ", "
                                       << "PDF set: " << steer<std::string>("set") << ", "
                                       << "member: " << steer<int>("member") << ".\n\t"
                                       << "x range: " << Limits{pdf_->xMin(), pdf_->xMax()} << ", "
                                       << "Q^2 range: " << Limits{pdf_->q2Min(), pdf_->q2Max()} << " GeV^2.\n\t"
                                       << "Interpolated from other flavours (" << pdf_->flavors()
                                       << "): " << from_remnant_ << ".";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("LHAPDF collinear photon flux");
        //desc.add<std::string>("set", "NNPDF31_nnlo_pdfas").setDescription("PDFset to use");
        //desc.add<std::string>("set", "LUXqed17_plus_PDF4LHC15_nnlo_100").setDescription("PDFset to use");
        //desc.add<std::string>("set", "LUXlep-NNPDF31_nlo_as_0118_luxqed").setDescription("PDFset to use");
        //desc.add<std::string>("set", "LUXqed_plus_PDF4LHC15_nnlo_100").setDescription("PDFset to use");
        desc.add<std::string>("set", "cteq66").setDescription("PDFset to use");
        desc.add<int>("member", 0).setDescription("PDF member");
        desc.addAs<int, pdgid_t>("pdgId", (pdgid_t)22).setDescription("parton PDG identifier");
        desc.add<bool>("fromRemnant", true)
            .setDescription("extrapolate distribution from sum imbalance of other contributions?");
        return desc;
      }

      double operator()(double x, double mx) const override {
        static const Limits x_valid{0., 1.};
        if (x == 0. || !x_valid.contains(x) || mx <= 0.)
          return 0.;
        const auto q2 = utils::q2(x, mp2_, mx * mx);
        if (!pdf_->inRangeXQ2(x, q2))
          return 0.;
        if (from_remnant_) {
          double xf = 0.;
          auto xf_flav = pdf_->xfxQ2(x, q2);
          for (const auto& flav : xf_flav)
            if (flav.first != (int)pdgid_)
              xf += flav.second;
          return prefactor_ * xf / x;
        }
        return prefactor_ * pdf_->xfxQ2((int)pdgid_, x, q2) / x;
      }

    private:
      std::unique_ptr<LHAPDF::PDF> pdf_;
      const pdgid_t pdgid_;
      const bool from_remnant_;
      static constexpr double prefactor_ = constants::ALPHA_EM * M_1_PI;
    };
  }  // namespace collflux
}  // namespace cepgen

typedef cepgen::collflux::LHAPDFCollinearFlux CF_LHAPDF;
REGISTER_COLLFLUX("LHAPDFCollinearFlux", CF_LHAPDF)
