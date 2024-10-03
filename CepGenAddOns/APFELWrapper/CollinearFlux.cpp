/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#include <APFEL/APFEL.h>

#include <cmath>

#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Message.h"

namespace cepgen::apfel {
  /// Generic partonic level perturbative structure functions built from an external PDFs grid
  class CollinearFlux final : public cepgen::CollinearFlux {
  public:
    /// Build a calculator from its Parameters object
    explicit CollinearFlux(const ParametersList& params)
        : cepgen::CollinearFlux(params), pdgid_(steer<pdgid_t>("partonPdgId")), q_range_(steer<Limits>("qrange")) {
      APFEL::SetPerturbativeOrder(steer<int>("perturbativeOrder"));
      const auto pdfset = steer<std::string>("set");
      if (!pdfset.empty())
        APFEL::SetPDFSet(pdfset);
      //APFEL::SetMaxFlavourPDFs(7);
      APFEL::SetFastEvolution(steer<bool>("fastEvolution"));
      APFEL::InitializeAPFEL();
      APFEL::EvolveAPFEL(q_range_.min(), q_range_.max());
      APFEL::CachePDFsAPFEL(q_range_.min());
      CG_INFO("apfel:CollinearFlux") << "Partonic collinear parton flux evaluator successfully built.\n"
                                     << " * APFEL version: " << APFEL::GetVersion() << "\n"
                                     << " * Parton PDG identifier: " << pdgid_
                                     << ", max flavours: " << APFEL::GetMaxFlavourPDFs() << "\n"
                                     << " * Q range: " << q_range_ << " (" << Limits{APFEL::GetMuF0(), APFEL::GetMuF()}
                                     << ") GeV\n"
                                     << " * perturbative order: " << APFEL::GetPerturbativeOrder() << ".";
    }

    cepgen::pdgid_t partonPdgId() const override { return pdgid_; }
    bool fragmenting() const override { return true; }
    double mass2() const override { return mp2_; }

    double fluxQ2(double x, double q2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      const auto q = std::sqrt(q2);
      if (!q_range_.contains(q))
        return 0.;
      return prefactor_ * APFEL::xPDFxQ(pdgid_, x, q);
    }

    static ParametersDescription description() {
      auto desc = cepgen::CollinearFlux::description();
      desc.setDescription("APFEL coll.flux");
      desc.add<std::string>("set", "").setDescription("LHAPDF set to use at the initial scale");
      desc.add<pdgid_t>("partonPdgId", PDG::photon).setDescription("parton PDG identifier");
      desc.add<Limits>("qrange", {1., 100.});
      desc.add<int>("perturbativeOrder", 2);
      desc.add<bool>("fastEvolution", false);
      return desc;
    }

  private:
    const pdgid_t pdgid_;
    const Limits q_range_;
  };
}  // namespace cepgen::apfel
using CollinearFluxAPFEL = cepgen::apfel::CollinearFlux;
REGISTER_COLLINEAR_FLUX("apfel", CollinearFluxAPFEL);
