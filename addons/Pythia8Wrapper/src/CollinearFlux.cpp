/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/PartonFluxes/CollinearFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"

using namespace std::string_literals;

namespace cepgen::pythia8 {
  class CollinearFlux final : public cepgen::CollinearFlux {
  public:
    explicit CollinearFlux(const ParametersList& params)
        : cepgen::CollinearFlux(params),
          type_(steer<std::string>("type")),
          parton_pdgid_(steer<pdgid_t>("partonPdgId")),
          beam_pdgid_(steer<pdgid_t>("beamPdgId")),
          mass2_(PDG::get().mass(beam_pdgid_)) {
      if (type_ == "Lepton") {
        const auto lepton_params = steer<ParametersList>("leptonParameters");
        info_.reset(new Pythia8::Info);
        if (const auto dil_sqrt_s = lepton_params.get<double>("sqrtS"); dil_sqrt_s > 0.)
          info_->setECM(dil_sqrt_s);
        else
          CG_WARNING("pythia8:CollinearFlux") << "Beam-beam centre-of-mass energy is required (through the 'sqrtS' "
                                                 "parameter) for the 'Lepton' collinear flux mode.";
        pdf_.reset(new Pythia8::Lepton(
            lepton_params.get<pdgid_t>("beamPdgId"), lepton_params.get<double>("Q2max"), info_.get()));
      } else if (type_ == "LHAGrid1")
        pdf_.reset(new Pythia8::LHAGrid1);
      else if (type_ == "MSTWpdf")
        pdf_.reset(new Pythia8::MSTWpdf);
      else if (type_ == "Proton2gammaDZ")
        pdf_.reset(new Pythia8::Proton2gammaDZ);
      else if (type_ == "Nucleus2gamma") {
        const auto hi_params = steer<ParametersList>("hiParameters");
        const auto hi = HeavyIon::fromPdgId(beam_pdgid_);
        const auto nucleon_mass = hi.mass(), b_min = hi_params.get<double>("bmin", hi.radius());
        pdf_.reset(new Pythia8::Nucleus2gamma(parton_pdgid_, b_min, nucleon_mass));
      } else if (type_ == "ProtonPoint")
        pdf_.reset(new Pythia8::ProtonPoint);
      else
        throw CG_FATAL("pythia8:CollinearFlux") << "Failed to initialise the Pythia 8 evaluator!\n"
                                                << "Parameters: " << params_;

      CG_INFO("pythia8:CollinearFlux") << "Pythia 8 '" << type_ << "' evaluator for collinear parton "
                                       << "(" << static_cast<PDG::Id>(beam_pdgid_) << " -> "
                                       << static_cast<PDG::Id>(parton_pdgid_) << ") flux initialised.";
    }

    static ParametersDescription description() {
      auto desc = cepgen::CollinearFlux::description();
      desc.setDescription("Pythia 8 coll.flux");
      desc.add("type", "Proton2gammaDZ"s)
          .allow("Lepton", "photon-from-lepton modelling")
          .allow("LHAGrid1", "LHAPDF grid modelling")
          .allow("MSTWpdf", "MSTW grid modelling")
          .allow("Proton2gammaDZ", "Drees-Zeppenfeld photon emission from proton")
          .allow("Nucleus2gamma", "photon-from-HI emission")
          .allow("ProtonPoint", "point-like photon emission from proton")
          .setDescription("type of PDF evaluator to use");
      desc.addAs<pdgid_t>("partonPdgId", PDG::photon).setDescription("parton PDG identifier");
      desc.addAs<pdgid_t>("beamPdgId", PDG::proton).setDescription("beam particle PDG identifier");
      auto lepton_desc = ParametersDescription();
      lepton_desc.add("sqrtS", -1.);
      lepton_desc.add("Q2max", 50.);
      desc.add("leptonParameters", lepton_desc);
      auto hi_desc = ParametersDescription();
      hi_desc.add("bmin", 0.).setDescription("minimum impact parameter for integration");
      desc.add("hiParameters", hi_desc);
      return desc;
    }

    pdgid_t partonPdgId() const override { return parton_pdgid_; }
    bool fragmenting() const override { return true; }
    double mass2() const override { return mass2_; }

    double fluxQ2(double x, double q2) const override {
      if (x == 0. || x < pdf_->getXmin())
        return 0.;
      return pdf_->xf(parton_pdgid_, x, q2);
    }

  private:
    std::unique_ptr<Pythia8::PDF> pdf_;
    std::unique_ptr<Pythia8::Info> info_;
    const std::string type_;
    const int parton_pdgid_, beam_pdgid_;
    const double mass2_;
  };
}  // namespace cepgen::pythia8
using PythiaCollinearFlux = cepgen::pythia8::CollinearFlux;
REGISTER_COLLINEAR_FLUX("pythia8", PythiaCollinearFlux);
