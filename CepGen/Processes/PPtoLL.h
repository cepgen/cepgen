#ifndef CepGen_Processes_PPtoLL_h
#define CepGen_Processes_PPtoLL_h

#include "GenericKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    /// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow \ell^+\ell^-\f$ process using \f$k_T\f$-factorization approach
    class PPtoLL : public GenericKTProcess
    {
      public:
        PPtoLL();
        ProcessPtr clone() const override { return ProcessPtr( new PPtoLL( *this ) ); }

      private:
        void preparePhaseSpace() override;
        /// \note IncQQbar in pptoll
        double computeKTFactorisedMatrixElement() override;
        void fillCentralParticlesKinematics() override;

        /// Rapidity range for the outgoing leptons
        Limits rap_limits_;
        /// Rapidity of the first outgoing lepton
        double y1_;
        /// Rapidity of the first outgoing lepton
        double y2_;

        /// Transverse momentum difference for the two outgoing leptons
        double pt_diff_;

        /// Azimuthal angle difference for the two outgoing leptons
        double phi_pt_diff_;

        /// First outgoing lepton's momentum
        Particle::Momentum Pl1_;
        /// Second outgoing lepton's momentum
        Particle::Momentum Pl2_;
    };
  }
}

#endif

