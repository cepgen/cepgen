/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2023  Laurent Forthomme
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
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// \brief CLAS parameterisation for nucleon data at \f$Q^2\f$ > 0.5 GeV\f$^2\f$ and \f$x_{\rm Bj}\f$ > 0.15
    /// \note This code was provided on 2016-04-13 by Silvano Simula and reflects the parameterisation used in \cite Osipenko:2003bu (CLAS) and described in \cite Ricco:1998yr.
    class CLAS : public Parameterisation {
    public:
      explicit CLAS(const ParametersList& params) : Parameterisation(params), mpi0_(PDG::get().mass(PDG::piZero)) {
        const auto& model = steer<std::string>("model");
        if (model == "proton")
          mod_params_ = Parameters::standard_proton();
        else if (model == "neutron")
          mod_params_ = Parameters::standard_neutron();
        else if (model == "deuteron")
          mod_params_ = Parameters::standard_deuteron();
        else
          throw CG_FATAL("CLAS") << "Invalid modelling selected: " << model << "!";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("CLAS (nucleon data, Q^2 > 0.5 GeV2 / xBj > 0.15)");
        desc.add<std::string>("model", "proton")
            .setDescription("Nucleon modelling ('proton', 'deuteron', or 'neutron' handled)");
        return desc;
      }

      /// List of steering parameters for a physics case
      struct Parameters {
        /// Standard parameterisation of a parton-from-neutron emission
        static Parameters standard_neutron();
        /// Standard parameterisation of a parton-from-proton emission
        static Parameters standard_proton();
        /// Standard parameterisation of a parton-from-deuteron emission
        static Parameters standard_deuteron();

        /// Physical properties associated to a resonance
        struct Resonance {
          double amplitude, mass, width;
          short angular_momentum;
        };

        enum { neutron = 0, proton = 1, deuteron = 2 } mode = proton;  ///< Nucleon type
        // SLAC fit parameters
        std::array<double, 7> c_slac = {};
        // CLAS parameterisation
        double alpha = 0., beta = 0., mu = 0., mup = 0.;
        std::array<double, 3> x = {};
        std::array<double, 4> b = {};
        std::vector<Resonance> resonances;
      };

      void eval() override {
        const double w2 = utils::mX2(args_.xbj, args_.q2, mp2_), w = sqrt(w2);
        if (w < mx_min_)
          return;
        const auto rb = resbkg(args_.q2, w);
        setF2(f2slac(args_.xbj, args_.q2) * (rb.first + rb.second));
      }

    private:
      /// \brief Method to evaluate the background/resonance terms of
      ///  the modulating function for the nucleon
      /// \note SLAC parameterisation
      std::pair<double, double> resbkg(double q2, double w) const;
      /// \brief Method to evaluate the deep inelastic structure function
      /// \f$F_{2}^{N}\f$ using the SLAC parameterisation
      /// \param[in] q2 squared four-momentum transfer in GeV\f$^2\f$
      /// \param[in] xbj Bjorken scaling variable
      /// \return \f$F_{2}^{N}\f$
      double f2slac(double xbj, double q2) const;
      Parameters mod_params_;
      static constexpr double COEFF = 6.08974;
      double mpi0_;  ///< Neutral pion mass
    };

    CLAS::Parameters CLAS::Parameters::standard_proton() {
      Parameters params;
      params.mode = Parameters::proton;
      // SLAC fit parameters
      params.c_slac = {{0.25615, 2.1785, 0.89784, -6.7162, 3.7557, 1.6421, 0.37636}};
      // CLAS parameterisation
      params.x = {{-0.599937, 4.76158, 0.411676}};
      params.b = {{0.755311, 3.35065, 3.51024, 1.74470}};
      params.alpha = -0.174985;
      params.beta = 0.00967019;
      params.mu = -0.0352567;
      params.mup = 3.51852;

      Parameters::Resonance r0{};
      r0.amplitude = 1.04;
      r0.mass = 1.22991;
      r0.width = 0.106254;
      r0.angular_momentum = 1;
      params.resonances.emplace_back(r0);

      Parameters::Resonance r1{};
      r1.amplitude = 0.481327;
      r1.mass = 1.51015;
      r1.width = 0.0816620;
      r1.angular_momentum = 2;
      params.resonances.emplace_back(r1);

      Parameters::Resonance r2{};
      r2.amplitude = 0.655872;
      r2.mass = 1.71762;
      r2.width = 0.125520;
      r2.angular_momentum = 3;
      params.resonances.emplace_back(r2);

      Parameters::Resonance r3{};
      r3.amplitude = 0.747338;
      r3.mass = 1.95381;
      r3.width = 0.198915;
      r3.angular_momentum = 2;
      params.resonances.emplace_back(r3);

      return params;
    }

    CLAS::Parameters CLAS::Parameters::standard_neutron() {
      Parameters params = standard_proton();
      params.mode = Parameters::neutron;
      params.c_slac = {{0.0640, 0.2250, 4.1060, -7.0790, 3.0550, 1.6421, 0.37636}};
      return params;
    }

    CLAS::Parameters CLAS::Parameters::standard_deuteron() {
      Parameters params = standard_proton();
      params.mode = Parameters::deuteron;
      params.c_slac = {{0.47709, 2.1602, 3.6274, -10.470, 4.9272, 1.5121, 0.35115}};
      params.x = {{-0.21262, 6.9690, 0.40314}};
      params.b = {{0.76111, 4.1470, 3.7119, 1.4218}};
      params.alpha = -0.24480;
      params.beta = 0.014503;

      params.resonances.clear();

      Parameters::Resonance r0{};
      r0.amplitude = 0.74847;
      r0.mass = 1.2400;
      r0.width = 0.12115;
      r0.angular_momentum = 1;
      params.resonances.emplace_back(r0);

      Parameters::Resonance r1{};
      r1.amplitude = 0.011500;
      r1.mass = 1.4772;
      r1.width = 0.0069580;
      r1.angular_momentum = 2;
      params.resonances.emplace_back(r1);

      Parameters::Resonance r2{};
      r2.amplitude = 0.12662;
      r2.mass = 1.5233;
      r2.width = 0.084095;
      r2.angular_momentum = 3;
      params.resonances.emplace_back(r2);

      Parameters::Resonance r3{};
      r3.amplitude = 0.747338;
      r3.mass = 1.95381;
      r3.width = 0.198915;
      r3.angular_momentum = 2;
      params.resonances.emplace_back(r3);

      return params;
    }

    double CLAS::f2slac(double xbj, double q2) const {
      if (xbj >= 1.)
        return 0.;

      const double xsxb = (q2 + mod_params_.c_slac[6]) / (q2 + mod_params_.c_slac[5] * xbj);
      const double xs = xbj * xsxb;

      double f2 = 0.;
      for (unsigned short i = 0; i < 5; ++i)
        f2 += mod_params_.c_slac[i] * pow(1. - xs, i);

      if (mod_params_.mode == Parameters::deuteron && xbj > 0.)
        f2 /= (1. - exp(-7.70 * (1. / xbj - 1. + mp2_ / q2)));

      return f2 * pow(1. - xs, 3) / xsxb;
    }

    std::pair<double, double> CLAS::resbkg(double q2, double w) const {
      const double mpi02 = mpi0_ * mpi0_;

      if (w < mx_min_)
        return std::make_pair(0., 0.);
      if (w > 4.)
        return std::make_pair(1., 0.);

      const double w2 = w * w;

      double qs = pow(w2 + mp2_ - mpi02, 2) - 4. * mp2_ * w2;
      if (qs <= 0.)
        return std::make_pair(1., 0.);
      qs = 0.5 * sqrt(qs) / w;

      const double omega = 0.5 * (w2 + q2 - mp2_) / mp_;
      const double xn = 0.5 * q2 / (mp_ * omega);

      const double bkg2 =
          (w > mod_params_.b[3]) ? exp(-mod_params_.b[2] * (w2 - mod_params_.b[3] * mod_params_.b[3])) : 1.;

      double f2bkg =
          (mod_params_.b[0]) * (1. - exp(-mod_params_.b[1] * (w - mx_min_))) + (1. - mod_params_.b[0]) * (1. - bkg2);
      f2bkg *= (1. + (1. - f2bkg) * (mod_params_.x[0] + mod_params_.x[1] * pow(xn - mod_params_.x[2], 2)));

      double etab = 1., etad = 1.;
      if (mod_params_.mode != Parameters::deuteron && q2 <= 2. && w <= 2.5) {
        etab = 1. - 2.5 * q2 * exp(-12.5 * q2 * q2 - 50. * (w - 1.325) * (w - 1.325));
        etad = 1. + 2.5 * q2 * exp(-12.5 * q2 * q2);
      }
      f2bkg *= etab;

      double f2resn = 0.;

      unsigned short i = 0;
      for (const auto& res : mod_params_.resonances) {
        const double ai = (i == 0)
                              ? etad * (res.amplitude + q2 * std::min(0., mod_params_.alpha + mod_params_.beta * q2))
                              : res.amplitude;
        const double dmi = (i == 2) ? res.mass * (1. + mod_params_.mu / (1. + mod_params_.mup * q2)) : res.mass;
        double qs0 = pow(dmi * dmi + mp2_ - mpi02, 2) - 4. * mp2_ * dmi * dmi;
        if (qs0 <= 0.)
          break;
        qs0 = 0.5 * sqrt(qs0) / dmi;
        int ji = 2 * res.angular_momentum;
        const double dg =
            0.5 * res.width * pow(qs / qs0, ji + 1) * (1. + pow(COEFF * qs0, ji)) / (1. + pow(COEFF * qs, ji));
        f2resn += ai * dg / ((w - dmi) * (w - dmi) + dg * dg);
        ++i;
      }
      f2resn *= 0.5 * (1. - mod_params_.b[0]) * bkg2 / mp_ * M_1_PI;

      return std::make_pair(f2bkg, f2resn);
    }
  }  // namespace strfun
}  // namespace cepgen
using cepgen::strfun::CLAS;
REGISTER_STRFUN(103, CLAS);
