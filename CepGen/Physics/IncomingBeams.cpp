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
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/Modes.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

namespace cepgen {
  IncomingBeams::IncomingBeams(const ParametersList& params) : SteeredObject(params) {
    ParametersList plist_pos, plist_neg;

    // positive-z incoming beam
    pdgid_t pos_pdg = steer<int>("beam1id");
    const int hi_Z1 = steer<int>("beam1Z");
    if (hi_Z1 != 0)
      pos_pdg = HeavyIon(steer<int>("beam1A"), (Element)hi_Z1);
    const auto& hi_beam1 = steer<std::vector<int> >("heavyIon1");
    if (hi_beam1.size() == 2)
      pos_pdg = HeavyIon{(unsigned short)hi_beam1.at(0), (Element)hi_beam1.at(1)};
    else if (!hi_beam1.empty())
      throw CG_FATAL("Kinematics") << "Invalid format for first incoming beam's HI specification!\n\t"
                                   << "A pair of (A,Z) is required, got " << hi_beam1 << ".";

    // negative-z incoming beam
    pdgid_t neg_pdg = steer<int>("beam2id");
    const int hi_A2 = steer<int>("beam2A");
    const int hi_Z2 = steer<int>("beam2Z");
    if (hi_Z2 != 0)
      neg_pdg = HeavyIon(hi_A2, (Element)hi_Z2);
    const auto& hi_beam2 = steer<std::vector<int> >("heavyIon2");
    if (hi_beam2.size() == 2)
      neg_pdg = HeavyIon{(unsigned short)hi_beam2.at(0), (Element)hi_beam2.at(1)};
    else if (!hi_beam2.empty())
      throw CG_FATAL("Kinematics") << "Invalid format for second incoming beam's HI specification!\n\t"
                                   << "A pair of (A,Z) is required, got " << hi_beam2 << ".";

    //----- combined two-beam system

    //--- beams PDG ids
    if (params_.has<std::vector<ParametersList> >("pdgIds")) {
      const auto& beams_pdg = steer<std::vector<ParametersList> >("pdgIds");
      if (beams_pdg.size() == 2) {
        pos_pdg = abs(beams_pdg.at(0).get<int>("pdgid"));
        neg_pdg = abs(beams_pdg.at(1).get<int>("pdgid"));
      } else if (!beams_pdg.empty())
        throw CG_FATAL("Kinematics") << "Invalid list of PDG ids retrieved for incoming beams:\n\t"
                                     << "2 PDG ids are expected, " << beams_pdg << " provided.";
    } else if (params_.has<std::vector<int> >("pdgIds")) {
      const auto& beams_pdg = steer<std::vector<int> >("pdgIds");
      if (beams_pdg.size() == 2) {
        pos_pdg = abs(beams_pdg.at(0));
        neg_pdg = abs(beams_pdg.at(1));
      } else if (!beams_pdg.empty())
        throw CG_FATAL("Kinematics") << "Invalid list of PDG ids retrieved for incoming beams:\n\t"
                                     << "2 PDG ids are expected, " << beams_pdg << " provided.";
    }

    plist_pos.set<int>("pdgId", pos_pdg);
    plist_neg.set<int>("pdgId", neg_pdg);

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
    //--- centre-of-mass energy
    const double sqrts = steer<double>("sqrtS"), cme = steer<double>("cmEnergy");
    if (sqrts > 0. || cme > 0.) {
      if (pos_pdg != neg_pdg)
        throw CG_FATAL("Kinematics") << "Trying to set √s with asymmetric beams"
                                     << " (" << pos_pdg << "/" << neg_pdg << ").\n"
                                     << "Please fill incoming beams objects manually!";
      const auto max_sqrts = std::max(sqrts, cme);
      p1z = +0.5 * max_sqrts;
      p2z = -0.5 * max_sqrts;
    }
    if (p1z * p2z < 0. && p1z < 0.)
      std::swap(p1z, p2z);
    plist_pos.set<double>("pz", +fabs(p1z));
    plist_neg.set<double>("pz", -fabs(p2z));

    //--- form factors
    const auto ff_mode = steer<std::string>("formFactors");
    if (!ff_mode.empty() || !form_factors_)
      form_factors_ =
          formfac::FormFactorsFactory::get().build(ff_mode.empty() ? formfac::gFFStandardDipoleHandler : ff_mode);

    const auto mode = steerAs<int, mode::Kinematics>("mode");
    if (mode != mode::Kinematics::invalid) {
      switch (mode) {
        case mode::Kinematics::ElasticElastic:
          plist_pos.setAs<int, Beam::Mode>("mode", Beam::Mode::ProtonElastic);
          plist_neg.setAs<int, Beam::Mode>("mode", Beam::Mode::ProtonElastic);
          break;
        case mode::Kinematics::ElasticInelastic:
          plist_pos.setAs<int, Beam::Mode>("mode", Beam::Mode::ProtonElastic);
          plist_neg.setAs<int, Beam::Mode>("mode", Beam::Mode::ProtonInelastic);
          break;
        case mode::Kinematics::InelasticElastic:
          plist_pos.setAs<int, Beam::Mode>("mode", Beam::Mode::ProtonInelastic);
          plist_neg.setAs<int, Beam::Mode>("mode", Beam::Mode::ProtonElastic);
          break;
        case mode::Kinematics::InelasticInelastic:
          plist_pos.setAs<int, Beam::Mode>("mode", Beam::Mode::ProtonInelastic);
          plist_neg.setAs<int, Beam::Mode>("mode", Beam::Mode::ProtonInelastic);
          break;
        default:
          throw CG_FATAL("Kinematics:IncomingBeams:mode") << "Unsupported kinematics mode: " << mode << "!";
      }
    }
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
        plist_pos.set<int>("ktFlux", kt_fluxes.at(0));
        plist_neg.set<int>("ktFlux", kt_fluxes.size() > 1 ? kt_fluxes.at(1) : kt_fluxes.at(0));
      }
    } else if (params_.has<int>("ktFluxes")) {
      const auto& ktfluxes = steerAs<int, Beam::KTFlux>("ktFluxes");
      if (ktfluxes != Beam::KTFlux::invalid) {
        plist_pos.set<int>("ktFlux", (int)ktfluxes);
        plist_neg.set<int>("ktFlux", (int)ktfluxes);
      }
    }
    CG_DEBUG("IncomingBeams") << "Will build the following incoming beams:\n* " << plist_pos << "\n* " << plist_neg
                              << ".";
    pos_beam_ = Beam(plist_pos);
    neg_beam_ = Beam(plist_neg);
  }

  const ParametersList& IncomingBeams::parameters() const {
    params_ = SteeredObject::parameters();
    if (str_fun_)
      params_.set<ParametersList>("structureFunctions", str_fun_->parameters());
    params_.setAs<int, mode::Kinematics>("mode", mode())
        .set<int>("beam1id", pos_beam_.pdgId())
        .set<double>("beam1pz", +pos_beam_.momentum().pz())
        .set<int>("beam2id", neg_beam_.pdgId())
        .set<double>("beam2pz", -neg_beam_.momentum().pz())
        .set<std::vector<int> >("ktFluxes", {(int)pos_beam_.ktFlux(), (int)neg_beam_.ktFlux()})
        .set<double>("sqrtS", sqrtS());
    const HeavyIon hi1(pos_beam_.pdgId()), hi2(neg_beam_.pdgId());
    if (hi1)
      params_.set<int>("beam1A", hi1.A).set<int>("beam1Z", (int)hi1.Z);
    if (hi2)
      params_.set<int>("beam2A", hi2.A).set<int>("beam2Z", (int)hi2.Z);
    return params_;
  }

  void IncomingBeams::setSqrtS(double sqs) {
    if (pos_beam_.pdgId() != neg_beam_.pdgId())
      throw CG_FATAL("Kinematics") << "Trying to set √s with asymmetric beams"
                                   << " (" << pos_beam_.pdgId() << "/" << neg_beam_.pdgId() << ").\n"
                                   << "Please fill incoming beams objects manually!";
    pos_beam_.setMomentum(Momentum::fromPxPyPzM(0.,
                                                0.,
                                                +0.5 * sqs,
                                                HeavyIon::isHI(pos_beam_.pdgId())
                                                    ? HeavyIon::mass(HeavyIon(pos_beam_.pdgId()))
                                                    : PDG::get().mass(pos_beam_.pdgId())));
    neg_beam_.setMomentum(Momentum::fromPxPyPzM(0.,
                                                0.,
                                                -0.5 * sqs,
                                                HeavyIon::isHI(neg_beam_.pdgId())
                                                    ? HeavyIon::mass(HeavyIon(neg_beam_.pdgId()))
                                                    : PDG::get().mass(neg_beam_.pdgId())));
  }

  double IncomingBeams::s() const { return (pos_beam_.momentum() + neg_beam_.momentum()).mass2(); }

  double IncomingBeams::sqrtS() const { return std::sqrt(s()); }

  mode::Kinematics IncomingBeams::mode() const {
    const auto& pos_mode = pos_beam_.parameters().getAs<int, Beam::Mode>("mode");
    const auto& neg_mode = neg_beam_.parameters().getAs<int, Beam::Mode>("mode");
    switch (pos_mode) {
      case Beam::Mode::PointLikeFermion:
      case Beam::Mode::ProtonElastic: {
        switch (neg_mode) {
          case Beam::Mode::ProtonElastic:
          case Beam::Mode::PointLikeFermion:
            return mode::Kinematics::ElasticElastic;
          default:
            return mode::Kinematics::ElasticInelastic;
        }
      }
      case Beam::Mode::ProtonInelastic: {
        switch (neg_mode) {
          case Beam::Mode::ProtonElastic:
          case Beam::Mode::PointLikeFermion:
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
    sf_desc.setName<int>(11);  // default is SY
    desc.add<ParametersDescription>("structureFunctions", sf_desc)
        .setDescription("Beam inelastic structure functions modelling");
    desc.add<int>("ktFluxes", -1).setDescription("kT-factorised fluxes modelling");
    return desc;
  }
}  // namespace cepgen
