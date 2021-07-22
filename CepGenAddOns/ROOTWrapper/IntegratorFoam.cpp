#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Parameters.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"

namespace cepgen {
  /// Foam general-purpose integration algorithm
  /// as developed by S. Jadach (Institute of Nuclear Physics, Krakow, PL)
  class IntegratorFoam final : public Integrator, public TFoamIntegrand {
  public:
    explicit IntegratorFoam(const ParametersList&);
    static std::string description() { return "FOAM general purpose MC integrator"; }

    void integrate(double&, double&) override;
    inline double uniform() const override { return rnd_->Rndm(); }

    /// Compute the weight for a given phase space point
    inline double Density(int ndim, double* x) override {
      if (!integrand_)
        throw CG_FATAL("FoamDensity") << "Integrand object not yet initialised!";
      return integrand_->eval(std::vector<double>(x, x + ndim));
    }

  private:
    std::unique_ptr<TFoam> foam_;
    std::unique_ptr<TRandom> rnd_;
  };

  IntegratorFoam::IntegratorFoam(const ParametersList& params) : Integrator(params), foam_(new TFoam("Foam")) {
    const auto& rnd_mode = params.get<std::string>("rngEngine", "MersenneTwister");
    if (rnd_mode == "Ranlux")
      rnd_.reset(new TRandom1);
    else if (rnd_mode == "generic")
      rnd_.reset(new TRandom2);
    else if (rnd_mode == "MersenneTwister")
      rnd_.reset(new TRandom3);
    else
      throw CG_FATAL("IntegratorFoam") << "Unrecognised random generator: \"" << rnd_mode << "\".";
    rnd_->SetSeed(seed_);

    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "FOAM integrator built\n\t"
                                 << "Version: " << foam_->GetVersion() << ".";
  }

  void IntegratorFoam::integrate(double& result, double& abs_error) {
    if (!initialised_) {
      foam_.reset(new TFoam("Foam"));
      foam_->SetPseRan(rnd_.get());
      foam_->SetnCells(params_.get<int>("nCells", 1000));
      foam_->SetnSampl(params_.get<int>("nSampl", 200));
      foam_->SetnBin(params_.get<int>("nBin", 8));
      foam_->SetEvPerBin(params_.get<int>("EvPerBin", 25));
      foam_->SetChat(std::max(verbosity_, 0));
      foam_->SetRho(this);
      foam_->SetkDim(integrand_->size());
      foam_->Initialize();
      initialised_ = true;
    }
    for (size_t i = 0; i < 100000; ++i)
      foam_->MakeEvent();
    //--- launch integration
    double norm, err;
    foam_->Finalize(norm, err);

    foam_->GetIntegMC(result, abs_error);
    result_ = result;
    err_result_ = abs_error;

    CG_DEBUG("IntegratorFoam").log([&](auto& log) {
      double eps = 5.e-4, avewt, wtmax, sigma;
      foam_->GetWtParams(eps, avewt, wtmax, sigma);
      const double ncalls = foam_->GetnCalls();
      const double effic = wtmax > 0 ? avewt / wtmax : 0.;
      log << "Result: " << result_ << " +- " << err_result_ << "\n\t"
          << "Relative error: " << err_result_ / result_ * 100. << "%\n\t"
          << "Dispersion/<wt>= " << sigma << ", <wt>= " << avewt << ", <wt>/wtmax= " << effic << ",\n\t"
          << " for epsilon = " << eps << "\n\t"
          << " nCalls (initialisation only)= " << ncalls << ".";
    });
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("Foam", IntegratorFoam)
