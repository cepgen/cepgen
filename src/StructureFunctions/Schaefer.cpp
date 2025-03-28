/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2025  Laurent Forthomme
 *                2017       Wolfgang Schaefer
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen::strfun {
  /// LUX-like hybrid modelling of \f$F_{2,L}\f$ structure functions
  class Schaefer final : public Parameterisation {
  public:
    /// User-steered Schäfer hybrid structure functions calculator
    explicit Schaefer(const ParametersList& params)
        : Parameterisation(params),
          q2_cut_(steer<double>("Q2cut")),
          w2_lim_(steer<std::vector<double> >("W2limits")),
          higher_twist_(steer<double>("higherTwist")),
          res_params_(steer<ParametersList>("resonancesSF")),
          pert_params_(steer<ParametersList>("perturbativeSF")),
          cont_params_(steer<ParametersList>("continuumSF")),
          resonances_model_(StructureFunctionsFactory::get().build(res_params_)),
          perturbative_model_(StructureFunctionsFactory::get().build(pert_params_)),
          continuum_model_(StructureFunctionsFactory::get().build(cont_params_)) {
      CG_DEBUG("LUXlike") << "LUXlike structure functions evaluator built with:\n"
                          << " * Q² cut:             " << q2_cut_ << " GeV²\n"
                          << " * W² ranges:          " << w2_lim_.at(0) << " GeV² / " << w2_lim_.at(1) << " GeV²\n"
                          << " *   resonances model: " << *resonances_model_ << "\n"
                          << " * perturbative model: " << *perturbative_model_ << "\n"
                          << " *    continuum model: " << *continuum_model_ << "\n"
                          << " * higher-twist corr:  " << higher_twist_ << ".";
      if (w2_lim_.size() < 2)
        throw CG_FATAL("LUXlike") << "Invalid number of transition regions for W^2. Should have two, got " << w2_lim_
                                  << ".";
      inv_omega_range_ = 1. / (w2_lim_.at(1) - w2_lim_.at(0));
      if (inv_omega_range_ <= 0.)
        throw CG_FATAL("LUXlike") << "Invalid W^2 transition regions definitions: " << w2_lim_.at(0) << " / "
                                  << w2_lim_.at(1) << " GeV^2!";
    }

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("LUXlike (hybrid)");
      desc.add("Q2cut", 9.);
      desc.add("W2limits", std::vector{3., 4.});
      desc.add("higherTwist", 5.5);
      desc.add("resonancesSF", StructureFunctionsFactory::get().describeParameters("ChristyBosted"));
      desc.add("perturbativeSF", StructureFunctionsFactory::get().describeParameters("MSTWGrid"));
      desc.add("continuumSF", StructureFunctionsFactory::get().describeParameters("GD11p"));
      return desc;
    }

    void eval() override {
      const double w2 = utils::mX2(args_.xbj, args_.q2, mp2_);

      if (args_.q2 < q2_cut_) {
        if (w2 < w2_lim_.at(0)) {
          setF2(resonances_model_->F2(args_.xbj, args_.q2));
          setFL(resonances_model_->FL(args_.xbj, args_.q2));
        } else if (w2 < w2_lim_.at(1)) {
          const double r = rho(w2);
          setF2(r * continuum_model_->F2(args_.xbj, args_.q2) + (1. - r) * resonances_model_->F2(args_.xbj, args_.q2));
          setFL(r * continuum_model_->FL(args_.xbj, args_.q2) + (1. - r) * resonances_model_->FL(args_.xbj, args_.q2));
        } else {
          setF2(continuum_model_->F2(args_.xbj, args_.q2));
          setFL(continuum_model_->FL(args_.xbj, args_.q2));
        }
      } else {
        if (w2 < w2_lim_.at(1)) {
          setF2(continuum_model_->F2(args_.xbj, args_.q2));
          setFL(continuum_model_->FL(args_.xbj, args_.q2));
        } else {
          setF2(perturbative_model_->F2(args_.xbj, args_.q2));
          setFL(perturbative_model_->FL(args_.xbj, args_.q2) * (1. + higher_twist_ / args_.q2));
        }
      }
    }

  private:
    inline double rho(double w2) const {
      const double omega = (w2 - w2_lim_.at(0)) * inv_omega_range_;
      const double omega2 = omega * omega;
      return 2. * omega2 - omega2 * omega2;
    }

    const double q2_cut_;  ///< Transition \f$Q^2\f$ before reaching the continuum/perturbative regions
    /// Transition \f$W^2\f$ between:
    /// - resonances and hybrid continuum/resonances low-\f$Q^2\f$ regions,
    /// - hybrid continuum/resonances and continuum low-\f$Q^2\f$ regions, or
    /// - continuum and perturbative high-\f$Q^2\f$ regions
    const std::vector<double> w2_lim_;
    const double higher_twist_;  ///< Value of the higher-twist correction
    const ParametersList res_params_, pert_params_, cont_params_;

    const std::unique_ptr<Parameterisation> resonances_model_;  ///< Res.-dominated region (low-\f$Q^2/W^2\f$) modelling
    const std::unique_ptr<Parameterisation> perturbative_model_;  ///< Perturbative region (high-\f$Q^2/W^2\f$) modelling
    const std::unique_ptr<Parameterisation> continuum_model_;     ///< Continuum regions modelling
    double inv_omega_range_{-1.};
  };
}  // namespace cepgen::strfun
using cepgen::strfun::Schaefer;
REGISTER_STRFUN("LUXLike", 301, Schaefer);
