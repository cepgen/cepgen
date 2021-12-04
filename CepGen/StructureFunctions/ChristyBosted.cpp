/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include <utility>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
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

      struct Parameters {
        static Parameters standard();

        struct Resonance {
          explicit Resonance();

          double kr() const;
          double ecmr(double m2) const;
          double kcmr() const { return ecmr(0.); }
          double pcmr(double m2) const { return sqrt(std::max(0., ecmr(m2) * ecmr(m2) - m2)); }

          /// Branching ratios container for resonance decay into single, double pion or eta states
          struct BR {
            BR() : singlepi(0.), doublepi(0.), eta(0.) {}
            BR(double singlepi, double doublepi, double eta) : singlepi(singlepi), doublepi(doublepi), eta(eta) {}
            bool valid() const { return (singlepi + doublepi + eta == 1.); }
            /// single pion branching ratio
            double singlepi;
            /// double pion branching ratio
            double doublepi;
            /// eta meson branching ratio
            double eta;
          } br;
          /// proton mass, in GeV/c2
          double mp;
          /// meson angular momentum
          double angular_momentum;
          /// damping parameter
          double x0;
          /// mass, in GeV/c2
          double mass;
          /// full width, in GeV
          double width;
          double A0_T;
          double A0_L;
          std::array<double, 5> fit_parameters;
        };
        struct Continuum {
          struct Direction {
            explicit Direction() : sig0(0.) {}
            Direction(double sig0, std::vector<double> params) : sig0(sig0), fit_parameters(std::move(params)) {}
            double sig0;
            std::vector<double> fit_parameters;
          };
          std::array<Direction, 2> transverse;
          std::array<Direction, 1> longitudinal;
        } continuum;
        double m0;
        std::vector<Resonance> resonances;
      } mod_params_;
      const double mpi_, mpi2_;
      const double meta_, meta2_;
    };

    ChristyBosted::ChristyBosted(const ParametersList& params)
        : Parameterisation(params),
          mpi_(PDG::get().mass(PDG::piZero)),
          mpi2_(mpi_ * mpi_),
          meta_(PDG::get().mass(PDG::eta)),
          meta2_(meta_ * meta_) {
      const auto& model = params.get<std::string>("model");
      if (model == "standard")
        mod_params_ = Parameters::standard();
      else
        throw CG_FATAL("ChristyBosted") << "Invalid modelling selected: " << model << "!";
    }

    double ChristyBosted::resmod507(const Polarisation& pol, double w2, double q2) const {
      const double w = sqrt(w2);

      const double xb = utils::xBj(q2, mp2_, w2);
      double m0 = 0., q20 = 0.;

      switch (pol) {
        case Polarisation::T: {  // transverse
          m0 = 0.125, q20 = 0.05;
        } break;
        case Polarisation::L: {  // longitudinal
          m0 = mod_params_.m0, q20 = 0.125;
        } break;
      }

      const double norm_q2 = 1. / 0.330 / 0.330;
      const double t = log(log((q2 + m0) * norm_q2) / log(m0 * norm_q2));

      //--- calculate kinematics needed for threshold relativistic B-W
      // equivalent photon energies
      const double k = 0.5 * (w2 - mp2_) / mp_;
      const double kcm = 0.5 * (w2 - mp2_) / w;

      const double epicm = 0.5 * (w2 + mpi2_ - mp2_) / w, ppicm = sqrt(std::max(0., epicm * epicm - mpi2_));
      const double epi2cm = 0.5 * (w2 + 4. * mpi2_ - mp2_) / w,
                   ppi2cm = sqrt(std::max(0., epi2cm * epi2cm - 4 * mpi2_));
      const double eetacm = 0.5 * (w2 + meta2_ - mp2_) / w, petacm = sqrt(std::max(0., eetacm * eetacm - meta2_));

      std::array<double, 7> width, height, pgam;
      for (unsigned short i = 0; i < 7; ++i) {
        const Parameters::Resonance& res = mod_params_.resonances[i];
        width[i] = res.width;

        //--- calculate partial widths
        //----- 1-pion decay mode
        const double x02 = res.x0 * res.x0;
        const double partial_width_singlepi =
            pow(ppicm / res.pcmr(mpi2_), 2. * res.angular_momentum + 1.) *
            pow((res.pcmr(mpi2_) * res.pcmr(mpi2_) + x02) / (ppicm * ppicm + x02), res.angular_momentum);
        //----- 2-pion decay mode
        const double partial_width_doublepi =
            pow(ppi2cm / res.pcmr(4. * mpi2_), 2. * (res.angular_momentum + 2.)) *
            pow((res.pcmr(4. * mpi2_) * res.pcmr(4. * mpi2_) + x02) / (ppi2cm * ppi2cm + x02),
                res.angular_momentum + 2) *
            w / res.mass;
        //----- eta decay mode (only for S11's)
        const double partial_width_eta =
            (res.br.eta == 0.)
                ? 0.
                : pow(petacm / res.pcmr(meta2_), 2. * res.angular_momentum + 1.) *
                      pow((res.pcmr(meta2_) * res.pcmr(meta2_) + x02) / (petacm * petacm + x02), res.angular_momentum);

        // virtual photon width
        pgam[i] = res.width * pow(kcm / res.kcmr(), 2) * (res.kcmr() * res.kcmr() + x02) / (kcm * kcm + x02);

        width[i] = (partial_width_singlepi * res.br.singlepi + partial_width_doublepi * res.br.doublepi +
                    partial_width_eta * res.br.eta) *
                   res.width;

        //--- resonance Q^2 dependence calculations

        switch (pol) {
          case Polarisation::T:
            height[i] = res.A0_T * (1. + res.fit_parameters[0] * q2 / (1. + res.fit_parameters[1] * q2)) /
                        pow(1. + q2 / 0.91, res.fit_parameters[2]);
            break;
          case Polarisation::L:
            height[i] = res.A0_L / (1. + res.fit_parameters[3] * q2) * q2 * exp(-q2 * res.fit_parameters[4]);
            break;
        }
        height[i] *= height[i];
      }

      //--- calculate Breit-Wigners for all resonances

      double sig_res = 0.;
      for (unsigned short i = 0; i < 7; ++i) {
        const Parameters::Resonance& res = mod_params_.resonances[i];
        const double mass2 = res.mass * res.mass, width2 = width[i] * width[i];
        const double sigr = height[i] * res.kr() / k * res.kcmr() / kcm / res.width *
                            (width[i] * pgam[i] / (pow(w2 - mass2, 2) + mass2 * width2));
        sig_res += sigr;
      }
      sig_res *= w;

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
              const double expo = mod_params_.continuum.transverse[i].fit_parameters[1] +
                                  mod_params_.continuum.transverse[i].fit_parameters[2] * q2 +
                                  mod_params_.continuum.transverse[i].fit_parameters[3] * q2 * q2;
              sig_nr += mod_params_.continuum.transverse[i].sig0 /
                        pow(q2 + mod_params_.continuum.transverse[i].fit_parameters[0], expo) * pow(wdif, i + 1.5);
            }
          }

          sig_nr *= xpr;
        } break;
        case Polarisation::L: {  // longitudinal
          for (unsigned short i = 0; i < 1; ++i) {
            const double expo = mod_params_.continuum.longitudinal[i].fit_parameters[0] +
                                mod_params_.continuum.longitudinal[i].fit_parameters[1];
            sig_nr += mod_params_.continuum.longitudinal[i].sig0 * pow(1. - xpr, expo) / (1. - xb) *
                      pow(q2 / (q2 + q20), mod_params_.continuum.longitudinal[i].fit_parameters[2]) / (q2 + q20) *
                      pow(xpr,
                          mod_params_.continuum.longitudinal[i].fit_parameters[3] +
                              mod_params_.continuum.longitudinal[i].fit_parameters[4] * t);
          }
        }
      }
      return sig_res + sig_nr;
    }

    ChristyBosted::Parameters ChristyBosted::Parameters::standard() {
      Parameters params;

      params.m0 = 4.2802;
      params.continuum.transverse = {{Continuum::Direction(246.06, {{0.067469, 1.3501, 0.12054, -0.0038495}}),
                                      Continuum::Direction(-89.360, {{0.20977, 1.5715, 0.090736, 0.010362}})}};
      params.continuum.longitudinal = {{Continuum::Direction(86.746, {{0., 4.0294, 3.1285, 0.33403, 4.9623}})}};

      {  //--- P33(1232)
        Resonance p33;
        p33.br = Resonance::BR(1., 0., 0.);
        p33.angular_momentum = 1.;
        //p33.x0 = 0.15;
        p33.x0 = 0.14462;
        p33.mass = 1.2298;
        p33.width = 0.13573;
        p33.fit_parameters = {{4.2291, 1.2598, 2.1242, 19.910, 0.22587}};
        p33.A0_T = 7.7805;
        p33.A0_L = 29.414;
        params.resonances.emplace_back(p33);
      }
      {  //--- S11(1535)
        Resonance s11_1535;
        s11_1535.br = Resonance::BR(0.45, 0.1, 0.45);
        s11_1535.angular_momentum = 0.;
        s11_1535.x0 = 0.215;
        s11_1535.mass = 1.5304;
        s11_1535.width = 0.220;
        s11_1535.fit_parameters = {{6823.2, 33521., 2.5686, 0., 0.}};
        s11_1535.A0_T = 6.3351;
        s11_1535.A0_L = 0.;
        params.resonances.emplace_back(s11_1535);
      }
      {  //--- D13(1520)
        Resonance d13;
        d13.br = Resonance::BR(0.65, 0.35, 0.);
        d13.angular_momentum = 2.;
        d13.x0 = 0.215;
        d13.mass = 1.5057;
        d13.width = 0.082956;
        d13.fit_parameters = {{21.240, 0.055746, 2.4886, 97.046, 0.31042}};
        d13.A0_T = 0.60347;
        d13.A0_L = 157.92;
        params.resonances.emplace_back(d13);
      }
      {  //--- F15(1680)
        Resonance f15;
        f15.br = Resonance::BR(0.65, 0.35, 0.);
        f15.angular_momentum = 3.;
        f15.x0 = 0.215;
        f15.mass = 1.6980;
        f15.width = 0.095782;
        f15.fit_parameters = {{-0.28789, 0.18607, 0.063534, 0.038200, 1.2182}};
        f15.A0_T = 2.3305;
        f15.A0_L = 4.2160;
        params.resonances.emplace_back(f15);
      }
      {  //--- S11(1650)
        Resonance s11_1650;
        s11_1650.br = Resonance::BR(0.4, 0.5, 0.1);
        s11_1650.angular_momentum = 0.;
        s11_1650.x0 = 0.215;
        s11_1650.mass = 1.6650;
        s11_1650.width = 0.10936;
        s11_1650.fit_parameters = {{-0.56175, 0.38964, 0.54883, 0.31393, 2.9997}};
        s11_1650.A0_T = 1.9790;
        s11_1650.A0_L = 13.764;
        params.resonances.emplace_back(s11_1650);
      }
      {  //--- P11(1440) roper
        Resonance p11;
        p11.br = Resonance::BR(0.65, 0.35, 0.);
        p11.angular_momentum = 1.;
        p11.x0 = 0.215;
        p11.mass = 1.4333;
        p11.width = 0.37944;
        p11.fit_parameters = {{46.213, 0.19221, 1.9141, 0.053743, 1.3091}};
        p11.A0_T = 0.022506;
        p11.A0_L = 5.5124;
        params.resonances.emplace_back(p11);
      }
      {  //--- F37(1950)
        Resonance f37;
        f37.br = Resonance::BR(0.5, 0.5, 0.);
        f37.angular_momentum = 3.;
        f37.x0 = 0.215;
        f37.mass = 1.9341;
        f37.width = 0.380;
        f37.fit_parameters = {{0., 0., 1., 1.8951, 0.51376}};
        f37.A0_T = 3.4187;
        f37.A0_L = 1.8951;
        params.resonances.emplace_back(f37);
      }

      return params;
    }

    ChristyBosted::Parameters::Resonance::Resonance()
        : mp(PDG::get().mass(PDG::proton)), angular_momentum(0.), x0(0.), mass(0.), width(0.), A0_T(0.), A0_L(0.) {}

    double ChristyBosted::Parameters::Resonance::kr() const { return 0.5 * (mass * mass - mp * mp) / mp; }

    double ChristyBosted::Parameters::Resonance::ecmr(double m2) const {
      if (mass == 0.)
        return 0.;
      return 0.5 * (mass * mass + m2 - mp * mp) / mass;
    }

    ChristyBosted& ChristyBosted::eval(double xbj, double q2) {
      const double w2 = utils::mX2(xbj, q2, mp2_);

      if (sqrt(w2) < mx_min_) {
        F2 = 0.;
        return *this;
      }

      //-----------------------------
      // modification of Christy-Bosted at large q2 as described in the LUXqed paper
      //-----------------------------
      const double q21 = 30., q20 = 8.;
      const double delq2 = q2 - q20;
      const double qq = q21 - q20;
      //------------------------------

      double q2_eff = q2, w2_eff = w2;
      if (q2 > q20) {
        q2_eff = q20 + delq2 / (1. + delq2 / qq);
        w2_eff = utils::mX2(xbj, q2_eff, mp2_);
      }
      const double sigT = resmod507(Polarisation::T, w2_eff, q2_eff);
      const double sigL = resmod507(Polarisation::L, w2_eff, q2_eff);

      F2 = prefactor_ * (1. - xbj) * q2_eff / (1 + tau(xbj, q2_eff)) * (sigT + sigL) / constants::GEVM2_TO_PB * 1.e6;
      if (q2 > q20)
        F2 *= q21 / (q21 + delq2);

      if (sigT != 0.)
        Parameterisation::computeFL(q2_eff, xbj, sigL / sigT);

      return *this;
    }

    ParametersDescription ChristyBosted::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Christy-Bosted F2/FL parameterisation of low-mass resonances");
      desc.add<std::string>("model", "standard").setDescription("Parameterisation set (only 'standard' handled)");
      CG_WARNING("") << desc.parameters();
      return desc;
    }
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(strfun::Type::ChristyBosted, ChristyBosted, strfun::ChristyBosted)
