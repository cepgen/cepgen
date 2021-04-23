#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Core/Exception.h"

#include <Math/IntegratorMultiDim.h>

namespace cepgen {
  /// FOAM general-purpose integration algorithm
  class IntegratorROOT : public Integrator {
  public:
    explicit IntegratorROOT(const ParametersList&);
    static std::string description() { return "ROOT general purpose MC integrator"; }

    void integrate(double&, double&) override;

  private:
    std::function<double(const double*)> func_;
    std::vector<double> min_, max_;
    /// integration type (adaptive, MC methods, etc..)
    const std::string type_;
    /// desired absolute Error
    const double absTol_;
    /// desired relative Error
    const double relTol_;
    /// maximum number of sub-intervals
    const unsigned int size_;

    std::unique_ptr<ROOT::Math::IntegratorMultiDim> integr_;
  };

  IntegratorROOT::IntegratorROOT(const ParametersList& params)
      : Integrator(params),
        func_([=](const double* x) -> double {
          return integrand_->eval(std::vector<double>(x, x + integrand_->size()));
        }),
        type_(params.get<std::string>("type", "default")),
        absTol_(params.get<double>("absTol", -1.)),
        relTol_(params.get<double>("relTol", -1.)),
        size_(params.get<int>("size", 0)) {
    ROOT::Math::IntegratorMultiDim::Type type;
    if (type_ == "default")
      type = ROOT::Math::IntegratorMultiDim::Type::kDEFAULT;
    else if (type_ == "adaptive")
      type = ROOT::Math::IntegratorMultiDim::Type::kADAPTIVE;
    else if (type_ == "plain")
      type = ROOT::Math::IntegratorMultiDim::Type::kPLAIN;
    else if (type_ == "miser")
      type = ROOT::Math::IntegratorMultiDim::Type::kMISER;
    else if (type_ == "vegas")
      type = ROOT::Math::IntegratorMultiDim::Type::kVEGAS;
    else
      throw CG_FATAL("Integrator:build") << "Invalid integrator type specified: \"" << type_ << "\" is not supported!";
    integr_.reset(new ROOT::Math::IntegratorMultiDim(type, absTol_, relTol_, size_));
    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "ROOT generic integrator built\n\t"
                                 << "Type: " << integr_->Name() << ",\n\t"
                                 << "Absolute tolerance: " << absTol_ << ",\n\t"
                                 << "Relative tolerance: " << relTol_ << ",\n\t"
                                 << "Number of sub-intervals: " << size_ << ".";
  }

  void IntegratorROOT::integrate(double& result, double& abserr) {
    if (!initialised_) {
      integr_->SetFunction(func_, integrand_->size());
      min_ = std::vector<double>(integrand_->size(), 0.);
      max_ = std::vector<double>(integrand_->size(), 1.);
      initialised_ = true;
    }

    result_ = result = integr_->Integral(min_.data(), max_.data());
    err_result_ = abserr = integr_->Error();
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("ROOT", IntegratorROOT)
