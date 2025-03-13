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

#include <Math/Integrator.h>
#include <Math/IntegratorMultiDim.h>

#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/Message.h"

namespace cepgen::root {
  /// ROOT general-purpose integration algorithm
  class Integrator final : public cepgen::Integrator {
  public:
    explicit Integrator(const ParametersList& params)
        : cepgen::Integrator(params),
          type_(steer<std::string>("type")),
          absolute_tolerance_(steer<double>("absTol")),
          relative_tolerance_(steer<double>("relTol")),
          size_(steer<int>("size")) {
      {
        auto type = ROOT::Math::IntegratorMultiDim::Type::kDEFAULT;
        if (type_ == "adaptive")
          type = ROOT::Math::IntegratorMultiDim::Type::kADAPTIVE;
        else if (type_ == "plain")
          type = ROOT::Math::IntegratorMultiDim::Type::kPLAIN;
        else if (type_ == "miser")
          type = ROOT::Math::IntegratorMultiDim::Type::kMISER;
        else if (type_ == "vegas")
          type = ROOT::Math::IntegratorMultiDim::Type::kVEGAS;
        integrator_.reset(new ROOT::Math::IntegratorMultiDim(type, absolute_tolerance_, relative_tolerance_, size_));
      }
      {
        auto type = ROOT::Math::IntegratorOneDim::Type::kDEFAULT;
        if (type_ == "gauss")
          type = ROOT::Math::IntegratorOneDim::Type::kGAUSS;
        else if (type_ == "legendre")
          type = ROOT::Math::IntegratorOneDim::Type::kLEGENDRE;
        else if (type_ == "adaptive")
          type = ROOT::Math::IntegratorOneDim::Type::kADAPTIVE;
        else if (type_ == "adaptiveSingular")
          type = ROOT::Math::IntegratorOneDim::Type::kADAPTIVESINGULAR;
        else if (type_ == "nonAdaptive")
          type = ROOT::Math::IntegratorOneDim::Type::kNONADAPTIVE;
        integrator_1d_.reset(new ROOT::Math::IntegratorOneDim(type, absolute_tolerance_, relative_tolerance_, size_));
      }
      //--- a bit of printout for debugging
      CG_DEBUG("Integrator:build") << "ROOT generic integrator built\n\t"
                                   << "N-dimensional type: " << integrator_->Name() << ",\n\t"
                                   << "1-dimensional type: " << integrator_1d_->Name() << ",\n\t"
                                   << "Absolute tolerance: " << absolute_tolerance_ << ",\n\t"
                                   << "Relative tolerance: " << relative_tolerance_ << ",\n\t"
                                   << "Number of sub-intervals: " << size_ << ".";
    }

    static ParametersDescription description() {
      auto desc = cepgen::Integrator::description();
      desc.setDescription("ROOT general purpose MC integrator");
      desc.add<std::string>("type", "default");
      desc.add<double>("absTol", -1.);
      desc.add<double>("relTol", -1.);
      desc.add<int>("size", 0);
      return desc;
    }

    void setLimits(const std::vector<Limits>& lims) override {
      cepgen::Integrator::setLimits(lims);
      x_low_.clear();
      x_high_.clear();
      for (const auto& lim : limits_) {
        x_low_.emplace_back(lim.min());
        x_high_.emplace_back(lim.max());
      }
    }
    Value integrate(Integrand& integrand) override {
      checkLimits(integrand);

      if (integrand.size() == 1) {
        auto funct = [&](double x) -> double { return integrand.eval(std::vector{x}); };
        integrator_1d_->SetFunction(funct);
        return Value{integrator_1d_->Integral(limits_.at(0).min(), limits_.at(0).max()), integrator_1d_->Error()};
      }
      auto funct = [&](const double* x) -> double { return integrand.eval(std::vector(x, x + integrand.size())); };
      integrator_->SetFunction(funct, integrand.size());
      return Value{integrator_->Integral(x_low_.data(), x_high_.data()), integrator_->Error()};
    }

  private:
    const std::string type_;           ///< integration type (adaptive, MC methods, etc...)
    const double absolute_tolerance_;  ///< desired absolute Error
    const double relative_tolerance_;  ///< desired relative Error
    const unsigned int size_;          ///< maximum number of sub-intervals

    std::vector<double> x_low_, x_high_;
    std::unique_ptr<ROOT::Math::IntegratorMultiDim> integrator_;
    std::unique_ptr<ROOT::Math::IntegratorOneDim> integrator_1d_;
  };
}  // namespace cepgen::root
using ROOTIntegrator = cepgen::root::Integrator;
REGISTER_INTEGRATOR("root", ROOTIntegrator);
