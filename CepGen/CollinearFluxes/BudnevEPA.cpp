/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
 *                2003  Maarten Boonekamp, Tibor Kucs
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

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  namespace collflux {
    class BudnevEPALepton final : public Parameterisation {
    public:
      explicit BudnevEPALepton(const ParametersList& params)
          : Parameterisation(params), ml2_(std::pow(PDG::get().mass(steer<int>("pdgId")), 2)) {
        CG_INFO("BudnevEPALepton") << "Budnev EPA for photon-from-lepton elastic limit (lepton: "
                                   << PDG::get().name(steer<int>("pdgId")) << ").\n\t "
                                   << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
      }

      bool fragmenting() const override final { return false; }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Budnev EPA for lepton");
        desc.add<int>("pdgId", 11).setDescription("lepton PDG id");
        return desc;
      }

      double operator()(double x, double = 0.) const override {
        if (x >= 1.)
          return 0.;
        double q2min = ml2_ * x * x / (1. - x);
        if (!q2_range_.contains(q2min))
          return 0.;
        return std::max(0.,
                        0.5 * constants::ALPHA_EM * M_1_PI *
                            (2. * ml2_ * x * (-1. / q2min + 1. / q2_range_.max()) +
                             (2. - 2. * x + x * x) / x * log(q2_range_.max() / q2min)));
      }

    private:
      const double ml2_;
    };

    class BudnevEPANucleon : public Parameterisation {
    public:
      explicit BudnevEPANucleon(const ParametersList& params)
          : Parameterisation(params), a_(steer<double>("a")), b_(steer<double>("b")), c_(steer<double>("c")) {}

      bool fragmenting() const override final { return false; }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb()).setDescription("type of heavy ion considered");
        desc.add<double>("a", 7.16);
        desc.add<double>("b", -3.96);
        desc.add<double>("c", 0.028);
        return desc;
      }

      double operator()(double x, double = 0.) const override {
        if (x >= 1.)
          return 0.;
        double qmi = mass2() * x * x / (1. - x);
        if (!q2_range_.contains(qmi))
          return 0.;
        return std::max(0.,
                        constants::ALPHA_EM * M_1_PI * (phi_f(x, q2_range_.max() / qscale_) - phi_f(x, qmi / qscale_)) *
                            (1 - x) / x);
      }

    protected:
      double phi_f(double x, double qq) const {
        const double qq1 = 1 + qq, y = x * x / (1 - x);
        double f = (1 + a_ * y) * (-log(qq1 / qq) + 1 / qq1 + 1 / (2 * qq1 * qq1) + 1 / (3 * qq1 * qq1 * qq1));
        f += (1 - b_) * y / (4 * qq * qq1 * qq1 * qq1);
        f += c_ * (1 + y / 4) *
             (log((qq1 - b_) / qq1) + b_ / qq1 + b_ * b_ / (2 * qq1 * qq1) + b_ * b_ * b_ / (3 * qq1 * qq1 * qq1));
        return f;
      }
      virtual double mass2() const = 0;

    private:
      const double a_, b_, c_;
    };

    struct BudnevEPAProton final : public BudnevEPANucleon {
      explicit BudnevEPAProton(const ParametersList& params) : BudnevEPANucleon(params) {
        CG_INFO("BudnevEPAProton") << "Budnev EPA for photon-from-proton elastic limit.\n\t"
                                   << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
      }

      double mass2() const override { return mp2_; }

      static ParametersDescription description() {
        auto desc = BudnevEPANucleon::description();
        desc.setDescription("Budnev EPA for proton");
        return desc;
      }
    };

    class BudnevEPAHI final : public BudnevEPANucleon {
    public:
      explicit BudnevEPAHI(const ParametersList& params)
          : BudnevEPANucleon(params), hi_(steerAs<pdgid_t, HeavyIon>("heavyIon")), mass2_(hi_.mass() * hi_.mass()) {
        CG_INFO("BudnevEPAHI") << "Budnev EPA for photon-from-heavy ion elastic limit (HI: " << hi_ << ").\n\t"
                               << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
        if (q2_range_.max() < q2max_min_) {
          q2_range_.max() = q2max_min_;
          CG_INFO("BudnevEPAHI") << "Increased maximal Q^2 value to " << q2_range_.max() << ".";
        }
      }

      double mass2() const override { return mass2_; }

      static ParametersDescription description() {
        auto desc = BudnevEPANucleon::description();
        desc.setDescription("Budnev EPA for heavy ion");
        desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::Pb()).setDescription("type of heavy ion considered");
        return desc;
      }

      double operator()(double x, double = 0.) const override {
        return (unsigned short)hi_.Z * BudnevEPANucleon::operator()(x);
      }

    private:
      const HeavyIon hi_;
      const double mass2_;
      const double q2max_min_{1.e4};
    };
  }  // namespace collflux
}  // namespace cepgen

typedef cepgen::collflux::BudnevEPALepton CF_Lepton;
typedef cepgen::collflux::BudnevEPAHI CF_HI;
typedef cepgen::collflux::BudnevEPAProton CF_Proton;
REGISTER_COLLFLUX("BudnevEPALepton", CF_Lepton)
REGISTER_COLLFLUX("BudnevEPAHI", CF_HI)
REGISTER_COLLFLUX("BudnevEPAProton", CF_Proton)
