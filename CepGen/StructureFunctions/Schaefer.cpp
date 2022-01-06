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

#include <memory>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// LUX-like hybrid modelling of \f$F_{2,L}\f$ structure functions
    class Schaefer final : public Parameterisation {
    public:
      /// User-steered Schäfer hybrid structure functions calculator
      explicit Schaefer(const ParametersList&);

      static ParametersDescription description();

      Schaefer& eval(double xbj, double q2) override;
      std::string describe() const override;

    private:
      double rho(double w2) const;
      void initialise();
      /// Transition \f$Q^2\f$ before reaching the continuum/perturbative regions
      double q2_cut_;
      /// Transition \f$W^2\f$ between:
      /// - resonances and hybrid continuum/resonances low-\f$Q^2\f$ regions,
      /// - hybrid continuum/resonances and continuum low-\f$Q^2\f$ regions, or
      /// - continuum and perturbative high-\f$Q^2\f$ regions
      std::vector<double> w2_lim_;
      /// Enable/disable the HT correction
      bool higher_twist_;
      /// Resonances-dominated region (low-\f$Q^2/W^2\f$) modelling
      std::shared_ptr<Parameterisation> resonances_model_;
      /// Perturbative region (high-\f$Q^2/W^2\f$) modelling
      std::shared_ptr<Parameterisation> perturbative_model_;
      /// Continuum regions modelling
      std::shared_ptr<Parameterisation> continuum_model_;
      bool initialised_{false};
      double inv_omega_range_{-1.};
    };

    Schaefer::Schaefer(const ParametersList& params)
        : Parameterisation(params),
          q2_cut_(params.get<double>("Q2cut")),
          w2_lim_(params.get<std::vector<double> >("W2limits")),
          higher_twist_(params.get<bool>("higherTwist")) {
      const auto& res_params = params_.get<ParametersList>("resonancesSF");
      const auto& pert_params = params_.get<ParametersList>("perturbativeSF");
      const auto& cont_params = params_.get<ParametersList>("continuumSF");
      CG_DEBUG("Schaefer") << "LUXlike structure functions built using:\n"
                           << " *)   resonances: " << res_params << ",\n"
                           << " *) perturbative: " << pert_params << ",\n"
                           << " *)    continuum: " << cont_params << ".";
      resonances_model_ = StructureFunctionsFactory::get().build(res_params);
      perturbative_model_ = StructureFunctionsFactory::get().build(pert_params);
      continuum_model_ = StructureFunctionsFactory::get().build(cont_params);
    }

    std::string Schaefer::describe() const {
      std::ostringstream os;
      os << "LUXlike{"
         << "r=" << resonances_model_->describe() << ","
         << "p=" << perturbative_model_->describe() << ","
         << "c=" << continuum_model_->describe();
      if (higher_twist_)
        os << ",HT";
      os << "}";
      return os.str();
    }

    void Schaefer::initialise() {
      CG_DEBUG("LUXlike") << "LUXlike structure functions evaluator successfully initialised.\n"
                          << " * Q² cut:             " << q2_cut_ << " GeV²\n"
                          << " * W² ranges:          " << w2_lim_.at(0) << " GeV² / " << w2_lim_.at(1) << " GeV²\n"
                          << " *   resonances model: " << *resonances_model_ << "\n"
                          << " * perturbative model: " << *perturbative_model_ << "\n"
                          << " *    continuum model: " << *continuum_model_ << "\n"
                          << " * higher-twist?       " << std::boolalpha << higher_twist_;
      inv_omega_range_ = 1. / (w2_lim_.at(1) - w2_lim_.at(0));
      initialised_ = true;
    }

    Schaefer& Schaefer::eval(double xbj, double q2) {
      if (!initialised_)
        initialise();

      const double w2 = utils::mX2(xbj, q2, mp2_);

      if (q2 < q2_cut_) {
        if (w2 < w2_lim_.at(0)) {
          auto sf = (*resonances_model_)(xbj, q2);
          sf.computeFL(xbj, q2);
          F2 = sf.F2;
          FL = sf.FL;
          return *this;
        } else if (w2 < w2_lim_.at(1)) {
          auto sf_r = (*resonances_model_)(xbj, q2);
          auto sf_c = (*continuum_model_)(xbj, q2);
          sf_r.computeFL(xbj, q2);
          sf_c.computeFL(xbj, q2);
          const double r = rho(w2);
          F2 = r * sf_c.F2 + (1. - r) * sf_r.F2;
          FL = r * sf_c.FL + (1. - r) * sf_r.FL;
          return *this;
        } else {
          auto sf = (*continuum_model_)(xbj, q2);
          sf.computeFL(xbj, q2);
          F2 = sf.F2;
          FL = sf.FL;
          return *this;
        }
      } else {
        if (w2 < w2_lim_.at(1)) {
          auto sf = (*continuum_model_)(xbj, q2);
          sf.computeFL(xbj, q2);
          F2 = sf.F2;
          FL = sf.FL;
          return *this;
        } else {
          auto sf_p = (*perturbative_model_)(xbj, q2);
          F2 = sf_p.F2;
          sf_p.computeFL(xbj, q2);
          FL = sf_p.FL;
          if (higher_twist_)
            F2 *= (1. + 5.5 / q2);
          return *this;
        }
      }
    }

    double Schaefer::rho(double w2) const {
      if (inv_omega_range_ <= 0.)
        throw CG_FATAL("LUXlike") << "Invalid W² limits: " << w2_lim_.at(0) << " / " << w2_lim_.at(1) << " GeV²!";
      const double omega = (w2 - w2_lim_.at(0)) * inv_omega_range_;
      const double omega2 = omega * omega;
      return 2. * omega2 - omega * omega;
    }

    ParametersDescription Schaefer::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("LUXlike structure functions");
      desc.add<double>("Q2cut", 9.);
      desc.add<std::vector<double> >("W2limits", {3., 4.});
      desc.add<bool>("higherTwist", true);
      desc.add<ParametersDescription>("resonancesSF",
                                      StructureFunctionsFactory::get().describeParameters((int)Type::ChristyBosted));
      desc.add<ParametersDescription>("perturbativeSF",
                                      StructureFunctionsFactory::get().describeParameters((int)Type::MSTWgrid));
      desc.add<ParametersDescription>("continuumSF",
                                      StructureFunctionsFactory::get().describeParameters((int)Type::GD11p));
      return desc;
    }
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(strfun::Type::Schaefer, Schaefer, strfun::Schaefer)
