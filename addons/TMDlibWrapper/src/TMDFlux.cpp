/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include <tmdlib/TMDlib.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/StreamCollector.h"

using namespace cepgen;
using namespace std::string_literals;

class TMDFlux final : public KTFlux {
public:
  explicit TMDFlux(const ParametersList& params) : KTFlux(params), parton_pdgid_(steer<int>("partonPdgId")) {
    tmd_.setVerbosity(steer<int>("verbosity"));
    if (const auto set = steer<std::string>("set"); !set.empty())
      if (const auto replica = steer<int>("replica"); replica >= 0)
        tmd_.TMDinit(set, replica);
      else
        tmd_.TMDinit(set);
    else
      throw CG_ERROR("TMDFlux") << "Failed to retrieve a set name.";
  }

  static ParametersDescription description() {
    auto desc = KTFlux::description();
    desc.setDescription("TMDlib kt-dependent flux");
    desc.add("set", "PB-NLO+QED-HERAI+II-set2"s).setDescription("dataset name");
    desc.add("verbosity", 99).setDescription("TMDlib evaluator verbosity");
    desc.add("replica", -1).setDescription("dataset replica");
    desc.addAs<int, pdgid_t>("partonPdgId", PDG::photon);
    return desc;
  }

  bool fragmenting() const override { return true; }
  double mass2() const override { return mp2_; }
  pdgid_t partonPdgId() const override { return parton_pdgid_; }

  double fluxQ2(double x, double kt2, double q2) const override {
    if (!utils::positive(x))
      return 0.;
    std::unordered_map<int, double> values_map;
    {
      std::string buf;
      auto sc = utils::StreamCollector(buf);
      tmd_.TMDpdf(x,
                  0.,  // xbar
                  std::sqrt(kt2),
                  std::sqrt(q2),  // evolution scale mu
                  values_map[PDG::up],
                  values_map[PDG::down],
                  values_map[4],  // sea
                  values_map[3],  // charm
                  values_map[5],  // bottom
                  values_map[PDG::gluon],
                  values_map[PDG::photon]);
    }
    if (values_map.count(parton_pdgid_) == 0)
      throw CG_ERROR("TMDFlux:fluxQ2") << "Parton id=" << static_cast<PDG::Id>(parton_pdgid_)
                                       << " is not handled by this TMD evaluator.";
    return values_map.at(parton_pdgid_);
  }

private:
  mutable TMDlib::TMD tmd_;
  const int parton_pdgid_;
};
REGISTER_KT_FLUX("tmd", 50, TMDFlux);
