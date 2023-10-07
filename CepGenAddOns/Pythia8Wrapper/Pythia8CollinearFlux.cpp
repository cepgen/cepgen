/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#include <Pythia8/PartonDistributions.h>

#include <numeric>

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  class Pythia8CollinearFlux : public CollinearFlux {
  public:
    explicit Pythia8CollinearFlux(const ParametersList& params)
        : CollinearFlux(params), type_(steer<std::string>("type")), pdgid_(steer<pdgid_t>("partonPdgId")) {
      if (type_ == "Lepton") {
        auto lepton_params = steer<ParametersList>("leptonParams");
        info_.reset(new Pythia8::Info);
        if (const auto dil_sqrt_s = lepton_params.get<double>("sqrtS"); dil_sqrt_s > 0.)
          info_->setECM(dil_sqrt_s);
        else
          CG_WARNING("Pythia8CollinearFlux") << "Beam-beam centre-of-mass energy is required (through the 'sqrtS' "
                                                "parameter) for the 'Lepton' collinear flux mode.";
        pdf_.reset(new Pythia8::Lepton(
            lepton_params.get<pdgid_t>("beamPdgId"), lepton_params.get<double>("Q2max"), info_.get()));
      } else if (type_ == "LHAGrid1")
        pdf_.reset(new Pythia8::LHAGrid1);
      else if (type_ == "MSTWpdf")
        pdf_.reset(new Pythia8::MSTWpdf);
      else if (type_ == "Proton2gammaDZ")
        pdf_.reset(new Pythia8::Proton2gammaDZ);
      else if (type_ == "ProtonPoint")
        pdf_.reset(new Pythia8::ProtonPoint);
      else
        throw CG_FATAL("Pythia8CollinearFlux") << "Failed to initialise the Pythia 8 evaluator!\n"
                                               << "Parameters: " << params_;

      CG_INFO("Pythia8CollinearFlux") << "Pythia 8 '" << type_ << "' evaluator for collinear parton "
                                      << "(" << (PDG::Id)pdgid_ << ") flux initialised.";
    }

    static ParametersDescription description() {
      auto desc = CollinearFlux::description();
      desc.setDescription("Pythia 8 coll.flux");
      desc.add<std::string>("type", "Proton2gammaDZ").setDescription("type of PDF evaluator to use");
      desc.add<pdgid_t>("partonPdgId", PDG::photon).setDescription("parton PDG identifier");
      auto lepton_desc = ParametersDescription();
      lepton_desc.add<pdgid_t>("beamPdgId", PDG::electron).setDescription("beam particle PDG identifier");
      lepton_desc.add<double>("sqrtS", -1.);
      lepton_desc.add<double>("Q2max", 50.);
      desc.add<ParametersDescription>("leptonParams", lepton_desc);
      return desc;
    }

    pdgid_t partonPdgId() const override final { return pdgid_; }
    double mass2() const override final { return mp2_; }

    double fluxQ2(double x, double q2) const override {
      if (x == 0. || x < pdf_->getXmin())
        return 0.;
      return pdf_->xf((int)pdgid_, x, q2);
    }

  private:
    std::unique_ptr<Pythia8::PDF> pdf_;
    std::unique_ptr<Pythia8::Info> info_;
    const std::string type_;
    const pdgid_t pdgid_;
  };
}  // namespace cepgen
REGISTER_COLLINEAR_FLUX("pythia8", Pythia8CollinearFlux);
