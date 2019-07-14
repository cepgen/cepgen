#ifndef CepGen_Processes_FortranKTProcess_h
#define CepGen_Processes_FortranKTProcess_h

#include "CepGen/Processes/GenericKTProcess.h"
#include <functional>

namespace cepgen
{
  namespace proc
  {
    /// Compute the matrix element for a generic \f$k_{\rm T}\f$-factorised process defined in a Fortran weighting function
    class FortranKTProcess : public GenericKTProcess
    {
      public:
        FortranKTProcess( const ParametersList& params, const char* name, const char* descr, std::function<double(void)> func );
        ProcessPtr clone( const ParametersList& /*params*/ ) const override { return ProcessPtr( new FortranKTProcess( *this ) ); }

        static std::unordered_map<std::string,ParametersList> kProcParameters;

      private:
        void preparePhaseSpace() override;
        double computeKTFactorisedMatrixElement() override;
        void fillCentralParticlesKinematics() override;

        std::function<double(void)> func_; ///< Function to be called for weight computation
        double y1_; ///< First outgoing particle rapidity
        double y2_; ///< Second outgoing particle rapidity
        double pt_diff_; ///< Transverse momentum balance between outgoing particles
        double phi_pt_diff_; ///< Azimutal angle difference between outgoing particles

        Particle::Momentum mom_ip1_; ///< First incoming beam momentum
        Particle::Momentum mom_ip2_; ///< Second incoming beam momentum
    };
  }
}

#endif
