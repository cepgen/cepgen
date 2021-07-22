#ifndef CepGen_Processes_FortranKTProcess_h
#define CepGen_Processes_FortranKTProcess_h

#include <functional>

#include "CepGen/Processes/KTProcess.h"

namespace cepgen {
  namespace proc {
    /// Compute the matrix element for a generic \f$k_{\rm T}\f$-factorised process defined in a Fortran weighting function
    class FortranKTProcess : public KTProcess {
    public:
      FortranKTProcess(const ParametersList& params, std::function<double(void)> func);
      ProcessPtr clone() const override { return ProcessPtr(new FortranKTProcess(*this)); }

      static ParametersList kProcParameters;

    private:
      void preparePhaseSpace() override;
      double computeKTFactorisedMatrixElement() override;
      void fillCentralParticlesKinematics() override;

      std::function<double(void)> func_;  ///< Function to be called for weight computation
      double y1_;                         ///< First outgoing particle rapidity
      double y2_;                         ///< Second outgoing particle rapidity
      double pt_diff_;                    ///< Transverse momentum balance between outgoing particles
      double phi_pt_diff_;                ///< Azimuthal angle difference between outgoing particles

      Momentum mom_ip1_;  ///< First incoming beam momentum
      Momentum mom_ip2_;  ///< Second incoming beam momentum
    };
  }  // namespace proc
}  // namespace cepgen

#endif
