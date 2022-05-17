/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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
#include "CepGen/Physics/Utils.h"

namespace cepgen {
  namespace collflux {
    class LHAPDFCollinearFlux : public Parameterisation {
    public:
      explicit LHAPDFCollinearFlux(const ParametersList& params)
          : Parameterisation(params), pdf_(LHAPDF::mkPDF(steer<std::string>("set"), steer<int>("member"))) {
        if (!pdf_)
          throw CG_FATAL("LHAPDFCollinearFlux") << "Failed to initialise the LHAPDF evaluator!\n\t"
                                                << "Parameters: " << params_;
        CG_INFO("LHAPDFCollinearFlux") << "LHAPDF evaluator for collinear photon flux initialised.\n\t"
                                       << "PDF set: " << steer<std::string>("set") << ", "
                                       << "member: " << steer<int>("member") << ".\n\t"
                                       << "x range: " << Limits{pdf_->xMin(), pdf_->xMax()} << ", "
                                       << "Q^2 range: " << Limits{pdf_->q2Min(), pdf_->q2Max()} << " GeV^2.";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("LHAPDF collinear photon flux");
        //desc.add<std::string>("set", "LUXqed_plus_PDF4LHC15_nnlo_100").setDescription("PDFset to use");
        desc.add<std::string>("set", "MMHT2015qed_nlo").setDescription("PDFset to use");
        desc.add<int>("member", 0).setDescription("PDF member");
        return desc;
      }

      double operator()(double x, double mx) const override {
        if (x <= 0. || mx <= 0.)
          return 0.;
        const auto q2 = utils::q2(x, mp2_, mx * mx);
        return pdf_->xfxQ2(22, x, q2) / x;
      }

    private:
      std::unique_ptr<LHAPDF::PDF> pdf_;
    };
  }  // namespace collflux
}  // namespace cepgen

REGISTER_COLLFLUX(LHAPDFCollinearFlux, collflux::LHAPDFCollinearFlux);
