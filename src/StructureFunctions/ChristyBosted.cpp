/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include <numeric>
#include <utility>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ResonanceObject.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen::strfun {
  /// \f$F_{2,L}\f$ parameterisation by Christy and Bosted \cite Bosted:2007xd
  class ChristyBosted final : public Parameterisation {
  public:
    explicit ChristyBosted(const ParametersList& params)
        : Parameterisation(params),
          m0_(steer<double>("m0")),
          mx2_min_(mx_min_ * mx_min_),
          q20_(steer<double>("q20")),
          q21_(steer<double>("q21")),
          mpi_(PDG::get().mass(PDG::piZero)),
          mpi2_(mpi_ * mpi_),
          meta_(PDG::get().mass(PDG::eta)),
          meta2_(meta_ * meta_) {
      for (const auto& res : steer<std::vector<ParametersList> >("resonances"))
        resonances_.emplace_back(res);
      const auto& cont = steer<std::vector<ParametersList> >("continuum");
      if (cont.size() != 3)
        throw CG_FATAL("ChristyBosted") << "Continuum should have its three directions parameterisation defined! Found "
                                        << cont.size() << ".";
      for (const auto& cnt : cont)
        continuum_.emplace_back(cnt);
    }

    static ParametersDescription description();

    void eval() override {
      const double w2 = utils::mX2(args_.xbj, args_.q2, mp2_);
      if (w2 < mx2_min_)
        return;

      //-----------------------------
      // modification of Christy-Bosted at large q2 as described in the LUXqed paper
      //-----------------------------
      const double delq2 = args_.q2 - q20_;
      //------------------------------

      double q2_eff = args_.q2, w2_eff = w2;
      if (args_.q2 > q20_) {
        q2_eff = q20_ + delq2 / (1. + delq2 / (q21_ - q20_));
        w2_eff = utils::mX2(args_.xbj, q2_eff, mp2_);
      }
      const double sigT = resmod507(transverse, w2_eff, q2_eff);
      const double sigL = resmod507(longitudinal, w2_eff, q2_eff);

      double f2 = prefactor_ * (1. - args_.xbj) * q2_eff / gamma2(args_.xbj, q2_eff) * (sigT + sigL) /
                  constants::GEVM2_TO_PB * 1.e6;
      if (args_.q2 > q20_)
        f2 *= q21_ / (q21_ + delq2);
      setF2(f2);

      if (sigT != 0.)
        Parameterisation::computeFL(args_.xbj, q2_eff, sigL / sigT);
    }
    //--- already computed internally during F2 computation
    ChristyBosted& computeFL(double, double) override { return *this; }
    ChristyBosted& computeFL(double, double, double) override { return *this; }

  private:
    enum Polarisation { longitudinal, transverse };
    static constexpr double prefactor_ = 0.25 * M_1_PI * M_1_PI / constants::ALPHA_EM;

    double resmod507(const Polarisation& pol, double w2, double q2) const;

    class Resonance final : public ResonanceObject, public SteeredObject<Resonance> {
    public:
      explicit Resonance(const ParametersList& params)
          : ResonanceObject(params),
            SteeredObject<Resonance>(params),
            a0t_(SteeredObject<Resonance>::steer<double>("A0T")),
            a0l_(SteeredObject<Resonance>::steer<double>("A0L")),
            fit_pars_(SteeredObject<Resonance>::steer<std::vector<double> >("fitParameters")) {
        if (fit_pars_.size() != 5)
          throw CG_FATAL("ChristyBosted:Resonance")
              << "Invalid fit parameters multiplicity! " << fit_pars_.size() << " != 5.";
      }

      static ParametersDescription description() {
        auto desc = ResonanceObject::description();
        desc.add("A0T", 0.);
        desc.add("A0L", 0.);
        desc.add("fitParameters", std::vector(5, 0.));
        return desc;
      }

      double sigma(const Polarisation& pol, const KinematicsBlock& kin) const {
        const auto partial_width = partialWidth(kin), mass2 = mass_ * mass_;
        return height(pol, kin.q2) * kr() / kin.k * kcmr() / kin.kcm / width_ *
               (partial_width * photonWidth(kin) / (pow(kin.w2 - mass2, 2) + mass2 * partial_width * partial_width));
      }

    private:
      /// resonance Q^2 dependence
      double height(const Polarisation& pol, double q2) const {
        switch (pol) {
          case transverse:
            return std::pow(a0t_ * (1. + fit_pars_.at(0) * q2 / (1. + fit_pars_.at(1) * q2)) /
                                std::pow(1. + q2 / 0.91, fit_pars_.at(2)),
                            2);
          case longitudinal:
            return std::pow(a0l_ / (1. + fit_pars_.at(3) * q2) * q2 * exp(-q2 * fit_pars_.at(4)), 2);
          default:
            throw CG_FATAL("ChristyBosted:Resonance") << "Invalid polarisation state: " << static_cast<int>(pol) << "!";
        }
      }

      const double a0t_;
      const double a0l_;
      const std::vector<double> fit_pars_;
    };
    /// Continuum parameterisation along one direction
    struct ContinuumDirection final : SteeredObject<ContinuumDirection> {
      explicit ContinuumDirection(const ParametersList& params)
          : SteeredObject(params), sig0(steer<double>("sig0")), fit_pars(steer<std::vector<double> >("fitParameters")) {
        if (fit_pars.size() != 4)
          throw CG_FATAL("ChristyBosted:ContinuumDirection")
              << "The templated fit for a continuum direction should have 4 parameters! Found " << fit_pars.size()
              << ".";
      }

      static ParametersDescription description() {
        auto desc = ParametersDescription();
        desc.setDescription("Set of parameters for one given direction");
        desc.add("sig0", 0.);
        desc.add("fitParameters", std::vector(4, 0.));
        return desc;
      }
      double sig0;
      std::vector<double> fit_pars;
    };
    const double m0_;
    const double mx2_min_;
    std::vector<Resonance> resonances_;          ///< Collection of resonance parameterisations
    std::vector<ContinuumDirection> continuum_;  ///< Three-dimensional parameterisation of the continuum
    const double q20_;
    const double q21_;
    const double mpi_, mpi2_;
    const double meta_, meta2_;
  };

  double ChristyBosted::resmod507(const Polarisation& pol, double w2, double q2) const {
    const double w = sqrt(w2), q20 = pol == transverse ? 0.05 : 0.125;

    //--- kinematics needed for threshold relativistic B-W
    const auto kin = Resonance::KinematicsBlock(w2, q2, mp2_, mpi2_, meta2_);

    //--- calculate Breit-Wigner for all resonances
    const auto sigma_resonant =
        w * std::accumulate(resonances_.begin(), resonances_.end(), 0., [&pol, &kin](auto sig, const auto& res) {
          return sig + res.sigma(pol, kin);
        });

    //--- non-resonant background calculation
    const auto xpr = 1. / (1. + (w2 - mx2_min_) / (q2 + q20));
    if (xpr > 1.)
      return 0.;

    double sigma_non_resonant = 0.;
    switch (pol) {
      case transverse: {
        if (const double mass_diff = w - mx_min_; mass_diff >= 0.) {
          for (unsigned short i = 0; i < 2; ++i) {
            const auto& dir = continuum_.at(i);  // direction for this continuum
            const double expo = dir.fit_pars.at(1) + dir.fit_pars.at(2) * q2 + dir.fit_pars.at(3) * q2 * q2;
            sigma_non_resonant += dir.sig0 / pow(q2 + dir.fit_pars.at(0), expo) * pow(mass_diff, i + 1.5);
          }
        }
        sigma_non_resonant *= xpr;
      } break;
      case longitudinal: {
        const auto& direction = continuum_.at(2);
        const double expo = direction.fit_pars.at(0);
        static constexpr double norm_q2 = 1. / 0.330 / 0.330;
        const double t = std::log(std::log((q2 + m0_) * norm_q2) / std::log(m0_ * norm_q2));
        sigma_non_resonant += direction.sig0 * pow(1. - xpr, expo) / (1. - utils::xBj(q2, mp2_, w2)) *
                              pow(q2 / (q2 + q20), direction.fit_pars.at(1)) / (q2 + q20) *
                              pow(xpr, direction.fit_pars.at(2) + direction.fit_pars.at(3) * t);
      }
    }
    return sigma_resonant + sigma_non_resonant;
  }

  ParametersDescription ChristyBosted::description() {
    auto desc = Parameterisation::description();
    desc.setDescription("Christy-Bosted (low-mass resonances)");
    desc.add("m0", 4.2802);
    desc.add("q20", 8.).setDescription("Q^2 scale for the modification of the parameterisation");
    desc.add("q21", 30.);
    desc.addParametersDescriptionVector(
        "continuum",
        ContinuumDirection::description(),
        {// transverse direction x
         ParametersList().set("sig0", 246.06).set("fitParameters", std::vector{0.067469, 1.3501, 0.12054, -0.0038495}),
         // transverse direction y
         ParametersList().set("sig0", -89.360).set("fitParameters", std::vector{0.20977, 1.5715, 0.090736, 0.010362}),
         // longitudinal direction z
         ParametersList().set("sig0", 86.746).set("fitParameters", std::vector{4.0294, 3.1285, 0.33403, 4.9623})});
    desc.addParametersDescriptionVector(
            "resonances",
            Resonance::description(),
            {//--- P33(1232)
             ParametersList()
                 .set("branchingRatios", ParametersList().set("singlePi", 1.).set("doublePi", 0.).set("eta", 0.))
                 .set("angularMomentum", 1)
                 .set("x0", 0.14462 /* 0.15 */)
                 .set("mass", 1.2298)
                 .set("width", 0.13573)
                 .set("fitParameters", std::vector{4.2291, 1.2598, 2.1242, 19.910, 0.22587})
                 .set("A0T", 7.7805)
                 .set("A0L", 29.414),
             //--- S11(1535)
             ParametersList()
                 .set("branchingRatios", ParametersList().set("singlePi", 0.45).set("doublePi", 0.1).set("eta", 0.45))
                 .set("angularMomentum", 0)
                 .set("x0", 0.215)
                 .set("mass", 1.5304)
                 .set("width", 0.220)
                 .set("fitParameters", std::vector{6823.2, 33521., 2.5686, 0., 0.})
                 .set("A0T", 6.3351)
                 .set("A0L", 0.),
             //--- D13(1520)
             ParametersList()
                 .set("branchingRatios", ParametersList().set("singlePi", 0.65).set("doublePi", 0.35).set("eta", 0.))
                 .set("angularMomentum", 2)
                 .set("x0", 0.215)
                 .set("mass", 1.5057)
                 .set("width", 0.082956)
                 .set("fitParameters", std::vector{21.240, 0.055746, 2.4886, 97.046, 0.31042})
                 .set("A0T", 0.60347)
                 .set("A0L", 157.92),
             //--- F15(1680)
             ParametersList()
                 .set("branchingRatios", ParametersList().set("singlePi", 0.65).set("doublePi", 0.35).set("eta", 0.))
                 .set("angularMomentum", 3)
                 .set("x0", 0.215)
                 .set("mass", 1.6980)
                 .set("width", 0.095782)
                 .set("fitParameters", std::vector{-0.28789, 0.18607, 0.063534, 0.038200, 1.2182})
                 .set("A0T", 2.3305)
                 .set("A0L", 4.2160),
             //--- S11(1650)
             ParametersList()
                 .set("branchingRatios", ParametersList().set("singlePi", 0.4).set("doublePi", 0.5).set("eta", 0.1))
                 .set("angularMomentum", 0)
                 .set("x0", 0.215)
                 .set("mass", 1.6650)
                 .set("width", 0.10936)
                 .set("fitParameters", std::vector{-0.56175, 0.38964, 0.54883, 0.31393, 2.9997})
                 .set("A0T", 1.9790)
                 .set("A0L", 13.764),
             //--- P11(1440) roper
             ParametersList()
                 .set("branchingRatios", ParametersList().set("singlePi", 0.65).set("doublePi", 0.35).set("eta", 0.))
                 .set("angularMomentum", 1)
                 .set("x0", 0.215)
                 .set("mass", 1.4333)
                 .set("width", 0.37944)
                 .set("fitParameters", std::vector{46.213, 0.19221, 1.9141, 0.053743, 1.3091})
                 .set("A0T", 0.022506)
                 .set("A0L", 5.5124),
             //--- F37(1950)
             ParametersList()
                 .set("branchingRatios", ParametersList().set("singlePi", 0.5).set("doublePi", 0.5).set("eta", 0.))
                 .set("angularMomentum", 3)
                 .set("x0", 0.215)
                 .set("mass", 1.9341)
                 .set("width", 0.380)
                 .set("fitParameters", std::vector{0., 0., 1., 1.8951, 0.51376})
                 .set("A0T", 3.4187)
                 .set("A0L", 1.8951)})
        .setDescription("collection of resonances modelled in this fit");

    return desc;
  }

}  // namespace cepgen::strfun
using cepgen::strfun::ChristyBosted;
REGISTER_STRFUN("ChristyBosted", 102, ChristyBosted);
