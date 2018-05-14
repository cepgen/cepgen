#ifndef CepGen_Processes_PPtoFF_h
#define CepGen_Processes_PPtoFF_h

#include "GenericKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    /// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow f\bar f\f$ process using \f$k_T\f$-factorization approach
    class PPtoFF : public GenericKTProcess
    {
      public:
        PPtoFF();
        ProcessPtr clone() const override { return ProcessPtr( new PPtoFF( *this ) ); }

      private:
        void preparePhaseSpace() override;
        /// \note IncQQbar in pptoll
        double computeKTFactorisedMatrixElement() override;
        void fillCentralParticlesKinematics() override;

        /// Rapidity range for the outgoing fermions
        double onShellME( double shat, double that, double uhat ) const;
        double offShellME( double shat, double that, double, double, double, double, const Particle::Momentum&, const Particle::Momentum&, const Particle::Momentum&, const Particle::Momentum& ) const;
        Limits rap_limits_;
        /// Rapidity of the first outgoing fermion
        double y1_;
        /// Rapidity of the first outgoing fermion
        double y2_;

        /// Transverse momentum difference for the two outgoing fermions
        double pt_diff_;

        /// Azimuthal angle difference for the two outgoing fermions
        double phi_pt_diff_;

        double mf_, mf2_, qf_;
        unsigned short colf_;
        /// First outgoing fermion's momentum
        Particle::Momentum p_f1_;
        /// Second outgoing fermion's momentum
        Particle::Momentum p_f2_;
    };
  }
}

#endif

