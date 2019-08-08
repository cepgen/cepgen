#ifndef CepGen_Processes_PPtoWW_h
#define CepGen_Processes_PPtoWW_h

#include "CepGen/Processes/Process2to4.h"

namespace cepgen
{
  namespace proc
  {
    /// \brief Compute the matrix element for a CE \f$\gamma\gamma\rightarrow W^+W^-\f$ process using \f$k_{\rm T}\f$-factorization approach
    /// \note The full theoretical description of this process definition may be found in \cite Luszczak:2018ntp.
    class PPtoWW : public Process2to4
    {
      public:
        PPtoWW( const ParametersList& params = ParametersList() );
        ProcessPtr clone( const ParametersList& params ) const override {
          return ProcessPtr( new PPtoWW( *this ) );
        }
        enum class Polarisation { full = 0, LL = 1, LT = 2, TL = 3, TT = 4 };

      private:
        static const double mw_, mw2_;

        void prepareKinematics() override;
        double computeCentralMatrixElement() const override;

        double amplitudeWW( double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4 ) const;
        double onShellME( double shat, double that, double uhat ) const;
        double offShellME( double shat, double that, double uhat, double phi_sum, double phi_diff ) const;

        const int method_;

        std::vector<short> pol_w1_, pol_w2_;
    };
  }
}

#endif

