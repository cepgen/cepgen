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

#include <array>
#include <numeric>
#include <utility>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ResonanceObject.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// \f$F_{2,L}\f$ parameterisation by Christy and Bosted \cite Bosted:2007xd
    class ChristyBosted final : public Parameterisation {
    public:
      explicit ChristyBosted(const ParametersList&);

      static ParametersDescription description();

      ChristyBosted& eval(double xbj, double q2) override;

      //--- already computed internally during F2 computation
      ChristyBosted& computeFL(double, double) override { return *this; }
      ChristyBosted& computeFL(double, double, double) override { return *this; }

    private:
      enum Polarisation { L, T };
      static constexpr double prefactor_ = 0.25 * M_1_PI * M_1_PI / constants::ALPHA_EM;

      double resmod507(const Polarisation& pol, double w2, double q2) const;

      class Resonance : public ResonanceObject, public SteeredObject<Resonance> {
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
          desc.add<double>("A0T", 0.);
          desc.add<double>("A0L", 0.);
          desc.add<std::vector<double> >("fitParameters", std::vector<double>(5, 0.));
          return desc;
        }

        double sigma(const Polarisation& pol, const KinematicsBlock& kin) const {
          const auto pwidth = partialWidth(kin), pwidth2 = pwidth * pwidth;
          const auto mass2 = mass_ * mass_;
          return height(pol, kin.q2) * kr() / kin.k * kcmr() / kin.kcm / width_ *
                 (pwidth * photonWidth(kin) / (pow(kin.w2 - mass2, 2) + mass2 * pwidth2));
        }

      private:
        /// resonance Q^2 dependence
        double height(const Polarisation& pol, double q2) const {
          switch (pol) {
            case Polarisation::T:
              return std::pow(a0t_ * (1. + fit_pars_.at(0) * q2 / (1. + fit_pars_.at(1) * q2)) /
                                  std::pow(1. + q2 / 0.91, fit_pars_.at(2)),
                              2);
            case Polarisation::L:
              return std::pow(a0l_ / (1. + fit_pars_.at(3) * q2) * q2 * exp(-q2 * fit_pars_.at(4)), 2);
            default:
              throw CG_FATAL("ChristyBosted:Resonance") << "Invalid polarisation state: " << (int)pol << "!";
          }
        }

        const double a0t_;
        const double a0l_;
        const std::vector<double> fit_pars_;
      };
      /// Continuum parameterisation along one direction
      struct ContinuumDirection : SteeredObject<ContinuumDirection> {
        explicit ContinuumDirection(const ParametersList& params)
            : SteeredObject(params),
              sig0(steer<double>("sig0")),
              fit_pars(steer<std::vector<double> >("fitParameters")) {
          if (fit_pars.size() != 4)
            throw CG_FATAL("ChristyBosted:ContinuumDirection")
                << "The templated fit for a continuum direction should have 4 parameters! Found " << fit_pars.size()
                << ".";
        }

        static ParametersDescription description() {
          auto desc = ParametersDescription();
          desc.setDescription("Set of parameters for one given direction");
          desc.add<double>("sig0", 0.);
          desc.add<std::vector<double> >("fitParameters", std::vector<double>(4, 0.));
          return desc;
        }
        double sig0;
        std::vector<double> fit_pars;
      };
      double m0_;
      /// Collection of resonance parameterisations
      std::vector<Resonance> resonances_;
      /// Three-dimensional parameterisation of the continuum
      std::vector<ContinuumDirection> continuum_;
      const double q20_;
      const double q21_;
      const double mpi_, mpi2_;
      const double meta_, meta2_;
    };

    ChristyBosted::ChristyBosted(const ParametersList& params)
        : Parameterisation(params),
          m0_(steer<double>("m0")),
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

    double ChristyBosted::resmod507(const Polarisation& pol, double w2, double q2) const {
      const double w = sqrt(w2);
      const double q20 = pol == Polarisation::T ? 0.05 : 0.125;

      //--- kinematics needed for threshold relativistic B-W
      const auto kin = Resonance::KinematicsBlock(w2, q2, mp2_, mpi2_, meta2_);

      //--- calculate Breit-Wigners for all resonances
      const auto sig_res =
          w * std::accumulate(resonances_.begin(), resonances_.end(), 0., [&pol, &kin](auto sig, const auto& res) {
            return sig + res.sigma(pol, kin);
          });

      //--- non-resonant background calculation
      const double xpr = 1. / (1. + (w2 - mx_min_ * mx_min_) / (q2 + q20));
      if (xpr > 1.)
        return 0.;

      double sig_nr = 0.;
      switch (pol) {
        case Polarisation::T: {  // transverse
          const double wdif = w - mx_min_;
          if (wdif >= 0.) {
            for (unsigned short i = 0; i < 2; ++i) {
              const auto& dir = continuum_.at(i);  // direction for this continuum
              const double expo = dir.fit_pars.at(1) + dir.fit_pars.at(2) * q2 + dir.fit_pars.at(3) * q2 * q2;
              sig_nr += dir.sig0 / pow(q2 + dir.fit_pars.at(0), expo) * pow(wdif, i + 1.5);
            }
          }

          sig_nr *= xpr;
        } break;
        case Polarisation::L: {  // longitudinal
          const auto& dir = continuum_.at(2);
          const double expo = dir.fit_pars.at(0);
          const double xb = utils::xBj(q2, mp2_, w2);
          const double norm_q2 = 1. / 0.330 / 0.330;
          const double t = log(log((q2 + m0_) * norm_q2) / log(m0_ * norm_q2));
          sig_nr += dir.sig0 * pow(1. - xpr, expo) / (1. - xb) * pow(q2 / (q2 + q20), dir.fit_pars.at(1)) / (q2 + q20) *
                    pow(xpr, dir.fit_pars.at(2) + dir.fit_pars.at(3) * t);
        }
      }
      return sig_res + sig_nr;
    }

    ChristyBosted& ChristyBosted::eval(double xbj, double q2) {
      const double w2 = utils::mX2(xbj, q2, mp2_);
      if (sqrt(w2) < mx_min_)
        return *this;

      //-----------------------------
      // modification of Christy-Bosted at large q2 as described in the LUXqed paper
      //-----------------------------
      const double delq2 = q2 - q20_;
      //------------------------------

      double q2_eff = q2, w2_eff = w2;
      if (q2 > q20_) {
        q2_eff = q20_ + delq2 / (1. + delq2 / (q21_ - q20_));
        w2_eff = utils::mX2(xbj, q2_eff, mp2_);
      }
      const double sigT = resmod507(Polarisation::T, w2_eff, q2_eff);
      const double sigL = resmod507(Polarisation::L, w2_eff, q2_eff);

      double f2 =
          prefactor_ * (1. - xbj) * q2_eff / gamma2(xbj, q2_eff) * (sigT + sigL) / constants::GEVM2_TO_PB * 1.e6;
      if (q2 > q20_)
        f2 *= q21_ / (q21_ + delq2);
      setF2(f2);

      if (sigT != 0.)
        Parameterisation::computeFL(q2_eff, xbj, sigL / sigT);

      return *this;
    }

    ParametersDescription ChristyBosted::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Christy-Bosted (low-mass resonances)");
      desc.add<double>("m0", 4.2802);
      desc.add<double>("q20", 8.).setDescription("Q^2 scale for the modification of the parameterisation");
      desc.add<double>("q21", 30.);
      desc.addParametersDescriptionVector(
          "continuum",
          ContinuumDirection::description(),
          {// transverse direction x
           ParametersList()
               .set<double>("sig0", 246.06)
               .set<std::vector<double> >("fitParameters", {0.067469, 1.3501, 0.12054, -0.0038495}),
           // transverse direction y
           ParametersList()
               .set<double>("sig0", -89.360)
               .set<std::vector<double> >("fitParameters", {0.20977, 1.5715, 0.090736, 0.010362}),
           // longitudinal direction z
           ParametersList()
               .set<double>("sig0", 86.746)
               .set<std::vector<double> >("fitParameters", {4.0294, 3.1285, 0.33403, 4.9623})});
      desc.addParametersDescriptionVector(
              "resonances",
              Resonance::description(),
              {//--- P33(1232)
               ParametersList()
                   .set<ParametersList>(
                       "branchingRatios",
                       ParametersList().set<double>("singlePi", 1.).set<double>("doublePi", 0.).set<double>("eta", 0.))
                   .set<int>("angularMomentum", 1)
                   .set<double>("x0", 0.14462 /* 0.15 */)
                   .set<double>("mass", 1.2298)
                   .set<double>("width", 0.13573)
                   .set<std::vector<double> >("fitParameters", {4.2291, 1.2598, 2.1242, 19.910, 0.22587})
                   .set<double>("A0T", 7.7805)
                   .set<double>("A0L", 29.414),
               //--- S11(1535)
               ParametersList()
                   .set<ParametersList>("branchingRatios",
                                        ParametersList()
                                            .set<double>("singlePi", 0.45)
                                            .set<double>("doublePi", 0.1)
                                            .set<double>("eta", 0.45))
                   .set<int>("angularMomentum", 0)
                   .set<double>("x0", 0.215)
                   .set<double>("mass", 1.5304)
                   .set<double>("width", 0.220)
                   .set<std::vector<double> >("fitParameters", {6823.2, 33521., 2.5686, 0., 0.})
                   .set<double>("A0T", 6.3351)
                   .set<double>("A0L", 0.),
               //--- D13(1520)
               ParametersList()
                   .set<ParametersList>("branchingRatios",
                                        ParametersList()
                                            .set<double>("singlePi", 0.65)
                                            .set<double>("doublePi", 0.35)
                                            .set<double>("eta", 0.))
                   .set<int>("angularMomentum", 2)
                   .set<double>("x0", 0.215)
                   .set<double>("mass", 1.5057)
                   .set<double>("width", 0.082956)
                   .set<std::vector<double> >("fitParameters", {21.240, 0.055746, 2.4886, 97.046, 0.31042})
                   .set<double>("A0T", 0.60347)
                   .set<double>("A0L", 157.92),
               //--- F15(1680)
               ParametersList()
                   .set<ParametersList>("branchingRatios",
                                        ParametersList()
                                            .set<double>("singlePi", 0.65)
                                            .set<double>("doublePi", 0.35)
                                            .set<double>("eta", 0.))
                   .set<int>("angularMomentum", 3)
                   .set<double>("x0", 0.215)
                   .set<double>("mass", 1.6980)
                   .set<double>("width", 0.095782)
                   .set<std::vector<double> >("fitParameters", {-0.28789, 0.18607, 0.063534, 0.038200, 1.2182})
                   .set<double>("A0T", 2.3305)
                   .set<double>("A0L", 4.2160),
               //--- S11(1650)
               ParametersList()
                   .set<ParametersList>(
                       "branchingRatios",
                       ParametersList().set<double>("singlePi", 0.4).set<double>("doublePi", 0.5).set<double>("eta", 0.1))
                   .set<int>("angularMomentum", 0)
                   .set<double>("x0", 0.215)
                   .set<double>("mass", 1.6650)
                   .set<double>("width", 0.10936)
                   .set<std::vector<double> >("fitParameters", {-0.56175, 0.38964, 0.54883, 0.31393, 2.9997})
                   .set<double>("A0T", 1.9790)
                   .set<double>("A0L", 13.764),
               //--- P11(1440) roper
               ParametersList()
                   .set<ParametersList>("branchingRatios",
                                        ParametersList()
                                            .set<double>("singlePi", 0.65)
                                            .set<double>("doublePi", 0.35)
                                            .set<double>("eta", 0.))
                   .set<int>("angularMomentum", 1)
                   .set<double>("x0", 0.215)
                   .set<double>("mass", 1.4333)
                   .set<double>("width", 0.37944)
                   .set<std::vector<double> >("fitParameters", {46.213, 0.19221, 1.9141, 0.053743, 1.3091})
                   .set<double>("A0T", 0.022506)
                   .set<double>("A0L", 5.5124),
               //--- F37(1950)
               ParametersList()
                   .set<ParametersList>(
                       "branchingRatios",
                       ParametersList().set<double>("singlePi", 0.5).set<double>("doublePi", 0.5).set<double>("eta", 0.))
                   .set<int>("angularMomentum", 3)
                   .set<double>("x0", 0.215)
                   .set<double>("mass", 1.9341)
                   .set<double>("width", 0.380)
                   .set<std::vector<double> >("fitParameters", {0., 0., 1., 1.8951, 0.51376})
                   .set<double>("A0T", 3.4187)
                   .set<double>("A0L", 1.8951)})
          .setDescription("collection of resonances modelled in this fit");

      return desc;
    }

  }  // namespace strfun
}  // namespace cepgen
typedef cepgen::strfun::ChristyBosted ChristyBosted;
REGISTER_STRFUN(102, ChristyBosted);
