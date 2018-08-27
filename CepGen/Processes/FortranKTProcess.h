#ifndef CepGen_Processes_FortranKTProcess_h
#define CepGen_Processes_FortranKTProcess_h

#include "GenericKTProcess.h"
#include <functional>

namespace CepGen
{
  namespace Process
  {
    /// Compute the matrix element for a generic \f$k_T\f$-factorised process defined in a Fortran subroutine
    class FortranKTProcess : public GenericKTProcess
    {
      public:
        FortranKTProcess( const ParametersList& params, const char* name, const char* descr, std::function<void(double&)> func );
        ProcessPtr clone() const override { return ProcessPtr( new FortranKTProcess( *this ) ); }

      private:
        void preparePhaseSpace() override;
        double computeKTFactorisedMatrixElement() override;
        void fillCentralParticlesKinematics() override;

        int pair_;
        int method_;
        /// Subroutine to be called for weight computation
        std::function<void(double&)> func_;
        double y1_, y2_, pt_diff_, phi_pt_diff_;

        Particle::Momentum mom_ip1_, mom_ip2_;
    };
  }
}

#endif

