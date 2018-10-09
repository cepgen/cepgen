#ifndef CepGen_Processes_PPtoWW_h
#define CepGen_Processes_PPtoWW_h

#include "CepGen/Processes/GenericKTProcess.h"
#include "CepGen/Core/ParametersList.h"

namespace CepGen
{
  namespace process
  {
    /// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_T\f$-factorization approach
    /// \note The full theoretical description of this process definition may be found in \cite Luszczak:2018ntp.
    class PPtoWW : public GenericKTProcess
    {
      public:
        PPtoWW( const ParametersList& params = ParametersList() );
        ProcessPtr clone( const ParametersList& params ) const override { return ProcessPtr( new PPtoWW( params ) ); }
        enum class Polarisation { full = 0, LL = 1, LT = 2, TL = 3, TT = 4 };

      private:
        static const double mw_, mw2_;

        void preparePhaseSpace() override;
        double computeKTFactorisedMatrixElement() override;
        void fillCentralParticlesKinematics() override;

        double amplitudeWW( double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4 );
        double onShellME( double shat, double that, double uhat );
        double offShellME( double shat, double that, double uhat, double phi_sum, double phi_diff );

        int method_;
        Polarisation pol_state_;
        std::vector<short> pol_w1_, pol_w2_;

        /// Rapidity range for the outgoing W bosons
        Limits rap_limits_;
        /// Rapidity of the first outgoing W boson
        double y1_;
        /// Rapidity of the first outgoing W boson
        double y2_;

        Limits ptdiff_limits_;
        /// Transverse momentum difference for the two outgoing W bosons
        double pt_diff_;

        Limits phi_pt_diff_limits_;
        /// Azimuthal angle difference for the two outgoing W bosons
        double phi_pt_diff_;

        /// First outgoing W boson's momentum
        Particle::Momentum p_w1_;
        /// Second outgoing W boson's momentum
        Particle::Momentum p_w2_;
    };
  }
}

#endif
