#include "CepGen/Core/Process2to4.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Processes/ProcessesHandler.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include <iomanip>

namespace cepgen
{
  namespace proc
  {
    /// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow f\bar f\f$ process using \f$k_{\rm T}\f$-factorization approach
    class PPtoFF : public Process2to4
    {
      public:
        PPtoFF( const ParametersList& params = ParametersList() );
        ProcessPtr clone( const ParametersList& params ) const override {
          return ProcessPtr( new PPtoFF( *this ) );
        }

      private:
        enum class ME { onShell = 0, offShell = 1 };
        void prepareKinematics() override;
        double computeCentralMatrixElement() const override;

        /// Rapidity range for the outgoing fermions
        double onShellME( double shat, double that, double uhat ) const;
        double offShellME() const;

        const ME method_;
        //--- parameters for off-shell matrix element
        unsigned short p_mat1_, p_mat2_;
        unsigned short p_term_ll_, p_term_lt_, p_term_tt1_, p_term_tt2_;

        double mf2_, qf_;
        double mA2_, mB2_;
        unsigned short colf_;
    };

    PPtoFF::PPtoFF( const ParametersList& params ) :
      Process2to4( params, "pptoff", "ɣɣ → f⁺f¯", { PDG::photon, PDG::photon }, params.get<ParticleProperties>( "pair", PDG::get()( PDG::muon ) ).pdgid ),
      method_ ( (ME)params.get<int>( "method", (int)ME::offShell ) ),
      p_mat1_( 0 ), p_mat2_( 0 ), p_term_ll_( 0 ), p_term_lt_( 0 ), p_term_tt1_( 0 ), p_term_tt2_( 0 ),
      mA2_( 0. ), mB2_( 0. )
    {
      if ( !cs_prop_.fermion || cs_prop_.charge == 0. )
        throw CG_FATAL( "PPtoFF:prepare" )
          << "Invalid fermion pair selected: " << cs_prop_.description
          << " (" << (int)cs_prop_.pdgid << ")!";

      mf2_ = cs_prop_.mass*cs_prop_.mass;
      qf_ = cs_prop_.charge/3.;
      colf_ = cs_prop_.colours;

      CG_DEBUG( "PPtoFF:prepare" )
        << "Produced particles: " << cs_prop_.description << " ("
        << "mass = " << cs_prop_.mass << " GeV, "
        << "charge = " << std::setprecision( 2 ) << qf_ << " e)\n\t"
        << "matrix element computation method: " << (int)method_ << ".";

      if ( method_ == ME::offShell ) { // off-shell matrix element
        const auto& ofp = params.get<ParametersList>( "offShellParameters" );
        p_mat1_ = ofp.get<int>( "mat1", 1 );
        p_mat2_ = ofp.get<int>( "mat2", 1 );
        p_term_ll_ = ofp.get<int>( "termLL", 1 );
        p_term_lt_ = ofp.get<int>( "termLT", 1 );
        p_term_tt1_ = ofp.get<int>( "termTT", 1 );
        p_term_tt2_ = ofp.get<int>( "termtt", 1 );
      }
    }

    void
    PPtoFF::prepareKinematics()
    {
      if ( !kin_.cuts.central.pt_diff.valid() )
        kin_.cuts.central.pt_diff = { 0., 50. }; // tighter cut for fermions
      mA2_ = (*event_)[Particle::IncomingBeam1][0].mass2();
      mB2_ = (*event_)[Particle::IncomingBeam2][0].mass2();
    }

    double
    PPtoFF::computeCentralMatrixElement() const
    {
      double amat2 = 0.;

      switch ( method_ ) {
        case ME::onShell: {
          //--- first compute Mendelstam variables
          const double shat = ( q1_+q2_ ).mass2();
          const double that1 = ( q1_-p_c1_ ).mass2(), that2 = ( q2_-p_c2_ ).mass2(), that = 0.5*( that1+that2 );
          const double uhat1 = ( q1_-p_c2_ ).mass2(), uhat2 = ( q2_-p_c1_ ).mass2(), uhat = 0.5*( uhat1+uhat2 );
          amat2 = onShellME( shat, that, uhat );
        } break;
        case ME::offShell: {
          amat2 = offShellME();
        } break;
      }

      const double g_em = 4.*M_PI*constants::ALPHA_EM*qf_*qf_;
      return amat2 * colf_ * ( g_em*g_em );
    }

    double
    PPtoFF::onShellME( double shat, double that, double uhat ) const
    {
      CG_DEBUG_LOOP( "PPtoFF:onShell" )
        << "shat: " << shat << ", that: " << that << ", uhat: " << uhat << ".";

      //=================================================================
      //     on-shell formula for M^2
      //=================================================================
      const double mf4 = mf2_*mf2_, mf8 = mf4*mf4;

      const double term1  =  6. *mf8,
                   term2  = -3. *mf4 *that*that,
                   term3  = -14.*mf4 *that*uhat,
                   term4  = -3. *mf4 *uhat*uhat,
                   term5  =      mf2_*that*that*that,
                   term6  =  7.* mf2_*that*that*uhat,
                   term7  =  7.* mf2_*that*uhat*uhat,
                   term8  =      mf2_*uhat*uhat*uhat,
                   term9  =          -that*that*that*uhat,
                   term10 =          -that*uhat*uhat*uhat;

      return -2.*( term1+term2+term3+term4+term5
                  +term6+term7+term8+term9+term10 )/( pow( ( mf2_-that )*( mf2_-uhat ), 2) );
    }

    double
    PPtoFF::offShellME() const
    {
      const double amt1 = std::hypot( p_c1_.pt(), cs_prop_.mass ), amt2 = std::hypot( p_c2_.pt(), cs_prop_.mass );
      const double alpha1 = amt1/sqs_*exp( y_c1_ ), beta1 = amt1/sqs_*exp( -y_c1_ );
      const double alpha2 = amt2/sqs_*exp( y_c2_ ), beta2 = amt2/sqs_*exp( -y_c2_ );
      const double x1 = alpha1+alpha2, x2 = beta1+beta2;
      const double z1p = alpha1/x1, z1m = alpha2/x1, z1 = z1p*z1m;
      const double z2p = beta1/x2, z2m = beta2/x2, z2 = z2p*z2m;

      CG_DEBUG_LOOP( "2to4:zeta" )
        << "z(1/2)p = " << z1p << " / " << z2p << "\n\t"
        << "z(1/2)m = " << z1m << " / " << z2m << ".";

      //--- photon virtualities
      const double t1abs = ( q1_.pt2() + x1*( MX_*MX_-mA2_ )+x1*x1*mA2_ )/( 1.-x1 );
      const double t2abs = ( q2_.pt2() + x2*( MY_*MY_-mB2_ )+x2*x2*mB2_ )/( 1.-x2 );

      CG_DEBUG_LOOP( "2to4:q2abs" ) << "Photon kinematics:\n\t"
        << "q1t2/q2t2 = " << q1_.pt2() << "/" << q2_.pt2() << "\n\t"
        << "t1/t2 = " << t1abs << "/" << t2abs << "\n\t"
        << "x1/x2 = " << x1 << "/" << x2 << ".";

      //--- epsilon-squared variable
      const double eps12 = mf2_+z1*t1abs, eps22 = mf2_+z2*t2abs;

      const Momentum ak1 = ( z1m*p_c1_-z1p*p_c2_ ), ak2 = ( z2m*p_c1_-z2p*p_c2_ );
      const Momentum ph_p1 = ak1+z1p*q2_, ph_m1 = ak1-z1m*q2_;
      const Momentum ph_p2 = ak2+z2p*q1_, ph_m2 = ak2-z2m*q1_;

      const double kp1 = 1./( ph_p1.pt2()+eps12 ), km1 = 1./( ph_m1.pt2()+eps12 );
      Momentum phi1 = kp1*ph_p1-km1*ph_m1;
      phi1.setPz( 0. );
      phi1.setEnergy( kp1-km1 );

      const double kp2 = 1./( ph_p2.pt2()+eps22 ), km2 = 1./( ph_m2.pt2()+eps22 );
      Momentum phi2 = kp2*ph_p2-km2*ph_m2;
      phi2.setPz( 0. );
      phi2.setEnergy( kp2-km2 );

      const double dot1 = phi1.threeProduct( q1_ )/qt1_, cross1 = phi1.crossProduct( q1_ )/qt1_;
      const double dot2 = phi2.threeProduct( q2_ )/qt2_, cross2 = phi2.crossProduct( q2_ )/qt2_;

      CG_DEBUG_LOOP( "PPtoFF:offShell" )
        << "q1 = " << q1_ << ", q1t = " << qt1_ << "\n\t"
        << "q2 = " << q2_ << ", q2t = " << qt2_ << "\n\t"
        << "phi1 = " << phi1 << "\n\t"
        << "phi2 = " << phi2 << "\n\t"
        << "(dot):   " << dot1 << " / " << dot2 << "\n\t"
        << "(cross): " << cross1 << " / " << cross2 << ".";

      const double aux2_1 = p_term_ll_ * ( mf2_ + 4.*z1*z1*t1abs ) * phi1.energy2()
                          + p_term_tt1_* ( ( z1p*z1p + z1m*z1m )*( dot1*dot1+cross1*cross1 ) )
                          + p_term_tt2_* ( cross1*cross1-dot1*dot1 )
                          - p_term_lt_ * 4.*z1*( z1p-z1m ) * phi1.energy() * qt1_*dot1;

      const double aux2_2 = p_term_ll_ * ( mf2_ + 4.*z2*z2*t2abs ) * phi2.energy2()
                          + p_term_tt1_* ( ( z2p*z2p + z2m*z2m )*( dot2*dot2+cross2*cross2 ) )
                          + p_term_tt2_* ( cross2*cross2-dot2*dot2 )
                          - p_term_lt_ * 4.*z2*( z2p-z2m ) * phi2.energy() * qt2_*dot2;

      //=================================================================
      //     convention of matrix element as in our kt-factorization
      //     for heavy flavours
      //=================================================================

      const double aux = 1./( q1_.pt2()*q2_.pt2() );
      const double amat2_1 = aux2_1*2.*z1*q1_.pt2()*aux,
                   amat2_2 = aux2_2*2.*z2*q2_.pt2()*aux;

      //=================================================================
      //     symmetrization
      //=================================================================

      const double amat2 = 0.5*( p_mat1_*amat2_1 + p_mat2_*amat2_2 );
      CG_DEBUG_LOOP( "PPtoFF:offShell" )
        << "aux2(1/2) = " << aux2_1 << " / " << aux2_2 << "\n\t"
        << "amat2(1/2), amat2 = " << amat2_1 << " / " << amat2_2 << " / " << amat2 << ".";

      return amat2*pow( x1*x2*s_, 2 );
    }
  }
}
// register process
REGISTER_PROCESS( "pptoff", PPtoFF )

