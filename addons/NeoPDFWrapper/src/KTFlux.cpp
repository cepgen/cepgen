/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#include <NeoPDF.hpp>
#include <cmath>

#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"

using namespace std::string_literals;

namespace cepgen::neopdf {
  class KTFlux : public cepgen::KTFlux {
  public:
    explicit KTFlux(const ParametersList& params)
        : cepgen::KTFlux(params),
          neopdf_(steer<std::string>("name"), steer<int>("member")),
          parton_pdg_id_(steer<int>("partonPdgId")) {}

    static ParametersDescription description() {
      auto desc = cepgen::KTFlux::description();
      desc.setDescription("NeoPDF kt-dependent flux");
      desc.add("name", "MAP22_grids_FF_Km_N3LL"s);
      desc.add("member", 0);
      desc.addAs<int>("partonPdgId", PDG::photon);
      return desc;
    }
    bool fragmenting() const final { return true; }
    double mass2() const override { return mp2_; }
    spdgid_t partonPdgId() const override { return parton_pdg_id_; }

    double fluxQ2(double x, double kt2, double q2) const override {
      if (x < neopdf_.x_min() || x > neopdf_.x_max())
        return 0.;
      if (q2 < neopdf_.q2_min() || q2 > neopdf_.q2_max())
        return 0.;
      return neopdf_.xfxQ2_ND(parton_pdg_id_, std::vector{std::sqrt(kt2), x, q2});  //FIXME
    }

  private:
    const ::neopdf::NeoPDF neopdf_;
    const spdgid_t parton_pdg_id_;
  };
}  // namespace cepgen::neopdf
using NeoPDFKTFlux = cepgen::neopdf::KTFlux;
REGISTER_KT_FLUX("NeoPDF", 200, NeoPDFKTFlux);
