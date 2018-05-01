#ifndef CepGen_Processes_PAtoLL_h
#define CepGen_Processes_PAtoLL_h

#include "GenericKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    /// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow \ell^+\ell^-\f$ process using \f$k_T\f$-factorization approach
    class PAtoLL : public GenericKTProcess
    {
      public:
        PAtoLL();
        ProcessPtr clone() const override { return ProcessPtr( new PAtoLL( *this ) ); }

      private:
        void preparePhaseSpace() override;
        /// \note IncQQbar in pptoll
        double computeKTFactorisedMatrixElement() override;
        void fillCentralParticlesKinematics() override;

        double y1_, y2_, pt_diff_, phi_pt_diff_;
    };
  }
}

#endif

