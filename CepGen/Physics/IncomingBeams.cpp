/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/Modes.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

namespace cepgen {
  IncomingBeams::IncomingBeams(const ParametersList& params) : SteeredObject(params) {
    // positive-z incoming beam
    positive().pdg = steer<int>("beam1id");
    const int hi_Z1 = steer<int>("beam1Z");
    if (hi_Z1 != 0)
      positive().pdg = HeavyIon(steer<int>("beam1A"), (Element)hi_Z1);
    const auto& hi_beam1 = steer<std::vector<int> >("heavyIon1");
    if (hi_beam1.size() == 2)
      positive().pdg = HeavyIon{(unsigned short)hi_beam1.at(0), (Element)hi_beam1.at(1)};
    else if (!hi_beam1.empty())
      throw CG_FATAL("Kinematics") << "Invalid format for first incoming beam's HI specification!\n\t"
                                   << "A pair of (A,Z) is required, got " << hi_beam1 << ".";

    // negative-z incoming beam
    negative().pdg = steer<int>("beam2id");
    const int hi_A2 = steer<int>("beam2A");
    const int hi_Z2 = steer<int>("beam2Z");
    if (hi_Z2 != 0)
      negative().pdg = HeavyIon(hi_A2, (Element)hi_Z2);
    const auto& hi_beam2 = steer<std::vector<int> >("heavyIon2");
    if (hi_beam2.size() == 2)
      negative().pdg = HeavyIon{(unsigned short)hi_beam2.at(0), (Element)hi_beam2.at(1)};
    else if (!hi_beam2.empty())
      throw CG_FATAL("Kinematics") << "Invalid format for second incoming beam's HI specification!\n\t"
                                   << "A pair of (A,Z) is required, got " << hi_beam2 << ".";

    //----- combined two-beam system

    //--- beams PDG ids
    if (params_.has<std::vector<ParametersList> >("pdgIds")) {
      const auto& beams_pdg = steer<std::vector<ParametersList> >("pdgIds");
      if (beams_pdg.size() == 2) {
        positive().pdg = abs(beams_pdg.at(0).get<int>("pdgid"));
        negative().pdg = abs(beams_pdg.at(1).get<int>("pdgid"));
      } else if (!beams_pdg.empty())
        throw CG_FATAL("Kinematics") << "Invalid list of PDG ids retrieved for incoming beams:\n\t"
                                     << "2 PDG ids are expected, " << beams_pdg << " provided.";
    } else if (params_.has<std::vector<int> >("pdgIds")) {
      const auto& beams_pdg = steer<std::vector<int> >("pdgIds");
      if (beams_pdg.size() == 2) {
        positive().pdg = abs(beams_pdg.at(0));
        negative().pdg = abs(beams_pdg.at(1));
      } else if (!beams_pdg.empty())
        throw CG_FATAL("Kinematics") << "Invalid list of PDG ids retrieved for incoming beams:\n\t"
                                     << "2 PDG ids are expected, " << beams_pdg << " provided.";
    }
    if (positive().pdg == PDG::electron)
      positive().mode = mode::Beam::PointLikeFermion;
    if (negative().pdg == PDG::electron)
      negative().mode = mode::Beam::PointLikeFermion;

    //--- beams longitudinal momentum
    double p1z = 0., p2z = 0;
    params_.fill<double>("beam1pz", p1z);
    params_.fill<double>("beam2pz", p2z);
    const auto& beams_pz = steer<std::vector<double> >("pz");
    if (beams_pz.size() == 2) {
      p1z = beams_pz.at(0);
      p2z = beams_pz.at(1);
    } else if (!beams_pz.empty())
      throw CG_FATAL("Kinematics") << "Invalid format for beams pz specification!\n\t"
                                   << "A vector of two pz's is required.";
    const HeavyIon hi1(positive().pdg), hi2(negative().pdg);
    const double m1 = hi1 ? HeavyIon::mass(hi1) : PDG::get().mass(positive().pdg);
    const double m2 = hi2 ? HeavyIon::mass(hi2) : PDG::get().mass(negative().pdg);

    if (p1z * p2z < 0. && p1z < 0.)
      std::swap(p1z, p2z);
    positive().momentum = Momentum::fromPxPyPzM(0., 0., +fabs(p1z), m1);
    negative().momentum = Momentum::fromPxPyPzM(0., 0., -fabs(p2z), m2);

    //--- centre-of-mass energy
    const double sqrts = steer<double>("sqrtS"), cme = steer<double>("cmEnergy");
    if (sqrts > 0. || cme > 0.)
      setSqrtS(std::max(sqrts, cme));
    //--- form factors
    const auto ff_mode = steer<std::string>("formFactors");
    if (!ff_mode.empty() || !form_factors_)
      form_factors_ =
          formfac::FormFactorsFactory::get().build(ff_mode.empty() ? formfac::gFFStandardDipoleHandler : ff_mode);

    const auto mode = steerAs<int, mode::Kinematics>("mode");
    if (mode != mode::Kinematics::invalid)
      setMode(mode);
    //--- structure functions
    auto strfun = steer<ParametersList>("structureFunctions");
    if (!strfun.empty() || !str_fun_) {
      CG_DEBUG("IncomingBeams") << "Structure functions modelling to be built: " << strfun << ".";
      str_fun_ = strfun::StructureFunctionsFactory::get().build(strfun);
    }
    //--- parton fluxes for kt-factorisation
    if (params_.has<std::vector<int> >("ktFluxes")) {
      auto kt_fluxes = steer<std::vector<int> >("ktFluxes");
      if (!kt_fluxes.empty()) {
        positive().kt_flux = (KTFlux)kt_fluxes.at(0);
        negative().kt_flux = (kt_fluxes.size() > 1) ? (KTFlux)kt_fluxes.at(1) : (KTFlux)kt_fluxes.at(0);
      }
    } else if (params_.has<int>("ktFluxes")) {
      const auto& ktfluxes = steerAs<int, KTFlux>("ktFluxes");
      if (ktfluxes != KTFlux::invalid)
        positive().kt_flux = negative().kt_flux = ktfluxes;
    }
  }

  const ParametersList& IncomingBeams::parameters() const {
    params_ = SteeredObject::parameters();
    if (str_fun_)
      params_.set<ParametersList>("structureFunctions", str_fun_->parameters());
    params_.set<int>("mode", (int)mode())
        .set<int>("beam1id", positive().pdg)
        .set<double>("beam1pz", +positive().momentum.pz())
        .set<int>("beam2id", negative().pdg)
        .set<double>("beam2pz", -negative().momentum.pz())
        .set<std::vector<int> >("ktFluxes", {(int)positive().kt_flux, (int)negative().kt_flux})
        .set<double>("sqrtS", sqrtS());
    const HeavyIon hi1(positive().pdg), hi2(negative().pdg);
    if (hi1)
      params_.set<int>("beam1A", hi1.A).set<int>("beam1Z", (int)hi1.Z);
    if (hi2)
      params_.set<int>("beam2A", hi2.A).set<int>("beam2Z", (int)hi2.Z);
    return params_;
  }

  void IncomingBeams::setSqrtS(double sqrts) {
    if (first.pdg != second.pdg)
      throw CG_FATAL("Kinematics") << "Trying to set √s with asymmetric beams"
                                   << " (" << first.pdg << "/" << second.pdg << ").\n"
                                   << "Please fill incoming beams objects manually!";
    positive().momentum = Momentum::fromPxPyPzM(0., 0., +0.5 * sqrts, PDG::get().mass(positive().pdg));
    negative().momentum = Momentum::fromPxPyPzM(0., 0., -0.5 * sqrts, PDG::get().mass(negative().pdg));
  }

  double IncomingBeams::s() const { return (positive().momentum + negative().momentum).mass2(); }

  double IncomingBeams::sqrtS() const { return std::sqrt(s()); }

  void IncomingBeams::setMode(const mode::Kinematics& mode) {
    switch (mode) {
      case mode::Kinematics::ElasticElastic:
        positive().mode = mode::Beam::ProtonElastic;
        negative().mode = mode::Beam::ProtonElastic;
        break;
      case mode::Kinematics::ElasticInelastic:
        positive().mode = mode::Beam::ProtonElastic;
        negative().mode = mode::Beam::ProtonInelastic;
        break;
      case mode::Kinematics::InelasticElastic:
        positive().mode = mode::Beam::ProtonInelastic;
        negative().mode = mode::Beam::ProtonElastic;
        break;
      case mode::Kinematics::InelasticInelastic:
        positive().mode = mode::Beam::ProtonInelastic;
        negative().mode = mode::Beam::ProtonInelastic;
        break;
      default:
        throw CG_FATAL("Kinematics:IncomingBeams:mode") << "Unsupported kinematics mode: " << mode << "!";
    }
  }

  mode::Kinematics IncomingBeams::mode() const {
    switch (positive().mode) {
      case mode::Beam::PointLikeFermion:
      case mode::Beam::ProtonElastic: {
        switch (negative().mode) {
          case mode::Beam::ProtonElastic:
          case mode::Beam::PointLikeFermion:
            return mode::Kinematics::ElasticElastic;
          default:
            return mode::Kinematics::ElasticInelastic;
        }
      }
      case mode::Beam::ProtonInelastic: {
        switch (negative().mode) {
          case mode::Beam::ProtonElastic:
          case mode::Beam::PointLikeFermion:
            return mode::Kinematics::InelasticElastic;
          default:
            return mode::Kinematics::InelasticInelastic;
        }
      }
      default:
        return mode::Kinematics::invalid;
    }
  }

  void IncomingBeams::setStructureFunctions(int sf_model, int sr_model) {
    const unsigned long kLHAPDFCodeDec = 10000000, kLHAPDFPartDec = 1000000;
    sf_model = (sf_model == 0 ? (int)strfun::Type::SuriYennie : sf_model);
    sr_model = (sr_model == 0 ? (int)sigrat::Type::SibirtsevBlunden : sr_model);
    auto sf_params = ParametersList().setName<int>(sf_model).set<int>("sigmaRatio", sr_model);
    if (sf_model / kLHAPDFCodeDec == 1) {  // SF from parton
      const unsigned long icode = sf_model % kLHAPDFCodeDec;
      sf_params.setName<int>((int)strfun::Type::Partonic)
          .set<int>("pdfId", icode % kLHAPDFPartDec)
          .set<int>("mode", icode / kLHAPDFPartDec);  // 0, 1, 2
    }
    CG_DEBUG("IncomingBeams:setStructureFunctions")
        << "Structure functions modelling to be built: " << sf_params << ".";
    str_fun_ = strfun::StructureFunctionsFactory::get().build(sf_params);
  }

  void IncomingBeams::setStructureFunctions(std::unique_ptr<strfun::Parameterisation> param) {
    str_fun_ = std::move(param);
  }

  ParametersDescription IncomingBeams::description() {
    auto desc = ParametersDescription();
    desc.add<int>("beam1id", 2212).setDescription("PDG id of the positive-z beam particle");
    desc.add<int>("beam1A", 1).setDescription("Atomic weight of the positive-z ion beam");
    desc.add<int>("beam1Z", 1).setDescription("Atomic number of the positive-z ion beam");
    desc.add<std::vector<int> >("heavyIon2", {}).setDescription("{A, Z} of the positive-z ion beam");
    desc.add<int>("beam2id", 2212).setDescription("PDG id of the negative-z beam particle");
    desc.add<int>("beam2A", 1).setDescription("Atomic weight of the negative-z ion beam");
    desc.add<int>("beam2Z", 1).setDescription("Atomic number of the negative-z ion beam");
    desc.add<std::vector<int> >("heavyIon2", {}).setDescription("{A, Z} of the negative-z ion beam");
    desc.add<std::vector<ParametersList> >("pdgIds", {}).setDescription("PDG description of incoming beam particles");
    desc.add<std::vector<int> >("pdgIds", {}).setDescription("PDG ids of incoming beam particles");
    desc.add<std::vector<double> >("pz", {}).setDescription("Beam momenta (in GeV/c)");
    desc.add<double>("sqrtS", -1.).setDescription("Two-beam centre of mass energy (in GeV)");
    desc.add<double>("cmEnergy", -1.).setDescription("Two-beam centre of mass energy (in GeV)");
    desc.add<std::string>("formFactors", "").setDescription("Beam form factors modelling");
    desc.add<int>("mode", (int)mode::Kinematics::invalid)
        .setDescription("Process kinematics mode (1 = elastic, (2-3) = single-dissociative, 4 = double-dissociative)");
    auto sf_desc = strfun::Parameterisation::description();
    sf_desc.add<int>(ParametersList::MODULE_NAME, 11);  // default is SY
    desc.add<ParametersDescription>("structureFunctions", sf_desc)
        .setDescription("Beam inelastic structure functions modelling");
    desc.add<int>("ktFluxes", -1).setDescription("kT-factorised fluxes modelling");
    return desc;
  }

  Beam::Beam() : pdg(PDG::proton), mode(mode::Beam::invalid), kt_flux(KTFlux::invalid) {}

  std::ostream& operator<<(std::ostream& os, const Beam& beam) {
    if ((HeavyIon)beam.pdg)
      os << (HeavyIon)beam.pdg;
    else
      os << (PDG::Id)beam.pdg;
    os << " (" << beam.momentum.pz() << " GeV/c), " << beam.mode;
    if (beam.kt_flux != KTFlux::invalid)
      os << " [unint.flux: " << beam.kt_flux << "]";
    return os;
  }
}  // namespace cepgen
