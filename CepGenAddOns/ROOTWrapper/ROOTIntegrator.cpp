/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

namespace cepgen {
  /// FOAM general-purpose integration algorithm
  class ROOTIntegrator final : public Integrator {
  public:
    explicit ROOTIntegrator(const ParametersList& params)
        : Integrator(params),
          type_(steer<std::string>("type")),
          absTol_(steer<double>("absTol")),
          relTol_(steer<double>("relTol")),
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
        integr_.reset(new ROOT::Math::IntegratorMultiDim(type, absTol_, relTol_, size_));
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
        integr_1d_.reset(new ROOT::Math::IntegratorOneDim(type, absTol_, relTol_, size_));
      }
      //--- a bit of printout for debugging
      CG_DEBUG("Integrator:build") << "ROOT generic integrator built\n\t"
                                   << "N-dimensional type: " << integr_->Name() << ",\n\t"
                                   << "1-dimensional type: " << integr_1d_->Name() << ",\n\t"
                                   << "Absolute tolerance: " << absTol_ << ",\n\t"
                                   << "Relative tolerance: " << relTol_ << ",\n\t"
                                   << "Number of sub-intervals: " << size_ << ".";
    }

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("ROOT general purpose MC integrator");
      desc.add<std::string>("type", "default");
      desc.add<double>("absTol", -1.);
      desc.add<double>("relTol", -1.);
      desc.add<int>("size", 0);
      return desc;
    }

    void setLimits(const std::vector<Limits>& lims) override {
      Integrator::setLimits(lims);
      xlow_.clear();
      xhigh_.clear();
      for (const auto& lim : limits_) {
        xlow_.emplace_back(lim.min());
        xhigh_.emplace_back(lim.max());
      }
    }
    Value integrate(Integrand& integrand) override {
      checkLimits(integrand);

      if (integrand.size() == 1) {
        auto funct = [&](double x) -> double { return integrand.eval(std::vector<double>{x}); };
        integr_1d_->SetFunction(funct);
        return Value{integr_1d_->Integral(limits_.at(0).min(), limits_.at(0).max()), integr_1d_->Error()};
      }
      auto funct = [&](const double* x) -> double {
        return integrand.eval(std::vector<double>(x, x + integrand.size()));
      };
      integr_->SetFunction(funct, integrand.size());
      return Value{integr_->Integral(xlow_.data(), xhigh_.data()), integr_->Error()};
    }

  private:
    const std::string type_;   ///< integration type (adaptive, MC methods, etc..)
    const double absTol_;      ///< desired absolute Error
    const double relTol_;      ///< desired relative Error
    const unsigned int size_;  ///< maximum number of sub-intervals

    std::vector<double> xlow_, xhigh_;
    std::unique_ptr<ROOT::Math::IntegratorMultiDim> integr_;
    std::unique_ptr<ROOT::Math::IntegratorOneDim> integr_1d_;
  };
}  // namespace cepgen
REGISTER_INTEGRATOR("root", ROOTIntegrator);
