#ifndef CepGen_Processes_FortranKTProcess_h
#define CepGen_Processes_FortranKTProcess_h

#include "GenericKTProcess.h"
#include <functional>

namespace CepGen
{
  namespace Process
  {
    /// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow \ell^+\ell^-\f$ process using \f$k_T\f$-factorization approach
    class FortranKTProcess : public GenericKTProcess
    {
      public:
        FortranKTProcess( const char* name, const char* descr, std::function<void(double&)> func );
        ProcessPtr clone() const override { return ProcessPtr( new FortranKTProcess( *this ) ); }

      private:
        void preparePhaseSpace() override;
        /// \note IncQQbar in pptoll
        double computeKTFactorisedMatrixElement() override;
        void fillCentralParticlesKinematics() override;

        std::function<void(double&)> func_;
        double y1_, y2_, pt_diff_, phi_pt_diff_;

        Particle::Momentum mom_ip1_, mom_ip2_;
    };
  }
}

#endif

