#include "CepGen/Physics/AlphaS.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

namespace {
  extern "C" {
  void initalphas_(int& iord, double& fr2, double& mur, double& asmur, double& mc, double& mb, double& mt);
  double alphas_(double& mur);
  }
}  // namespace

namespace cepgen {
  class AlphaSPEGASUS : public AlphaS {
  public:
    explicit AlphaSPEGASUS(const ParametersList& params)
        : AlphaS(params),
          iord_(params.get<int>("iord", 2)),
          fr2_(params.get<double>("fr2", 1.)),
          mur_(params.get<double>("mur", 1.)),
          asmur_(params.get<double>("asmur", 0.68183)) {
      double mc = PDG::get().mass(4), mb = PDG::get().mass(5), mt = PDG::get().mass(6);

      initalphas_(iord_, fr2_, mur_, asmur_, mc, mb, mt);
      CG_INFO("AlphaSPEGASUS:init") << "PEGASUS alpha(S) evolution algorithm initialised with parameters:\n\t"
                                    << "order: " << iord_ << ", fr2: " << fr2_ << ", "
                                    << "mur: " << mur_ << ", asmur: " << asmur_ << "\n\t"
                                    << "quark masses (GeV): charm: " << mc << ", bottom: " << mb << ", top: " << mt
                                    << ".";
    }
    static std::string description() { return "PEGASUS alphaS evolution algorithm"; }

    double operator()(double q) const override { return alphas_(q); }

  private:
    int iord_;
    double fr2_;
    double mur_;
    double asmur_;
  };
}  // namespace cepgen

REGISTER_ALPHAS_MODULE("pegasus", AlphaSPEGASUS)
