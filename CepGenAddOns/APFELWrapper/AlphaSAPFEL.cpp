#include "CepGen/Physics/AlphaS.h"
#include "CepGen/Core/Exception.h"

#include "APFEL/APFEL.h"

namespace cepgen {
  class AlphaSAPFEL : public AlphaS {
  public:
    explicit AlphaSAPFEL(const ParametersList& params)
        : AlphaS(params),
          order_(params.get<int>("order", 2)),
          q0_(params.get<double>("q0", 1.)),
          qmax_(params.get<double>("qmax", 10000.)) {
      APFEL::SetPerturbativeOrder(order_);
      APFEL::InitializeAPFEL();
      APFEL::EvolveAPFEL(q0_, qmax_);
    }
    static std::string description() { return "APFEL alphaS evolution algorithm"; }

    double operator()(double q) const override {
      if (q < q0_ || q > qmax_)
        CG_WARNING("AlphaSAPFEL:get") << "q = " << q << " outside the evolution range"
                                      << " [" << q0_ << ":" << qmax_ << "].";
      return APFEL::AlphaQCD(q);
    }

  private:
    int order_;
    double q0_, qmax_;
  };
}  // namespace cepgen

REGISTER_ALPHAS_MODULE("apfel", AlphaSAPFEL)
