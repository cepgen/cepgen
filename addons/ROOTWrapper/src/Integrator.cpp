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

#include "CepGen/Integration/BaseIntegrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/BaseIntegratorFactory.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace std::string_literals;

namespace cepgen::root {
  /// ROOT general-purpose integration algorithm
  class Integrator final : public cepgen::BaseIntegrator {
  public:
    explicit Integrator(const ParametersList& params)
        : cepgen::BaseIntegrator(params),
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
        integrator_1d_.reset(new ROOT::Math::IntegratorOneDim(
            type, absolute_tolerance_, relative_tolerance_, size_, steer<int>("rule")));
      }
      CG_DEBUG("Integrator:build") << "ROOT generic integrator built\n\t"
                                   << "N-dimensional type: " << integrator_->Name() << ",\n\t"
                                   << "1-dimensional type: " << integrator_1d_->Name() << ",\n\t"
                                   << "Absolute tolerance: " << absolute_tolerance_ << ",\n\t"
                                   << "Relative tolerance: " << relative_tolerance_ << ",\n\t"
                                   << "Number of sub-intervals: " << size_ << ".";
    }

    static ParametersDescription description() {
      auto desc = cepgen::BaseIntegrator::description();
      desc.setDescription("ROOT general purpose MC integrator");
      desc.add("type", "default"s).setDescription("type of integration");
      desc.add("absTol", -1.).setDescription("desired absolute error limit");
      desc.add("relTol", -1.).setDescription("desired relative error limit");
      desc.add("size", 0).setDescription("maximum number of sub-intervals to build");
      desc.add("rule", 0).setDescription("Gauss-Kronrod integration rule (only for GSL kADAPTIVE type)");
      return desc;
    }

    Value run(Integrand& integrand, const std::vector<Limits>& range) const override {
      if (integrand.size() == 1) {
        auto funct = [&](double x) -> double { return integrand.eval(std::vector{x}); };
        integrator_1d_->SetFunction(funct);
        return Value{integrator_1d_->Integral(range.at(0).min(), range.at(0).max()), integrator_1d_->Error()};
      }
      std::vector<double> x_low, x_high;
      for (const auto& dim_range : range) {
        x_low.emplace_back(dim_range.min());
        x_high.emplace_back(dim_range.max());
      }
      auto funct = [&](const double* x) -> double { return integrand.eval(std::vector(x, x + integrand.size())); };
      integrator_->SetFunction(funct, integrand.size());
      return Value{integrator_->Integral(x_low.data(), x_high.data()), integrator_->Error()};
    }

  private:
    const std::string type_;           ///< integration type (adaptive, MC methods, etc...)
    const double absolute_tolerance_;  ///< desired absolute Error
    const double relative_tolerance_;  ///< desired relative Error
    const unsigned int size_;          ///< maximum number of sub-intervals

    std::unique_ptr<ROOT::Math::IntegratorMultiDim> integrator_;
    std::unique_ptr<ROOT::Math::IntegratorOneDim> integrator_1d_;
  };
}  // namespace cepgen::root
using ROOTIntegrator = cepgen::root::Integrator;
REGISTER_BASE_INTEGRATOR("root", ROOTIntegrator);
