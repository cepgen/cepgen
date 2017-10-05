#include "ChristyBosted.h"

namespace CepGen
{
  namespace SF
  {
    double resmod507( int sf, double w2, double q2 )
    {
      const ChristyBostedParameterisation params = ChristyBostedParameterisation::standard();
      //struct tst2 { std::array<double,7> sigrsv; double sig_nrsv, sig_mec; };

      const double mp = Constants::mp, mp2 = mp*mp,
                   mpi = Constants::mpi, mpi2 = mpi*mpi,
                   meta = Particle::massFromPDGId( Particle::Eta ), meta2 = meta*meta;
      const double w = sqrt( w2 );

      const double mon = 1./( 1.+q2/0.71 ),
                   dip = mon*mon; // dipole parameterization

      const double xb = q2/( q2+w2-mp2 );
      double m0 = 0., q20 = 0., xpr = 0.;

      if ( sf == 1 ) { // transverse
        m0 = 0.125, q20 = 0.05;
        xpr = 1./( 1.+( w2-pow( mp+mpi, 2 ) )/( q2+q20 ) );
      }
      else if ( sf == 2 ) {
        m0 = params.m0, q20 = 0.125;
        xpr = 1./( 1.+( w2-pow( mp+mpi+mpi, 2 ) )/( q2+q20 ) );
      }

      const double norm_q2 = 1./0.330/0.330;
      const double t = log( log( ( q2+m0 )*norm_q2 )/log( m0*norm_q2 ) );

      //--- calculate kinematics needed for threshold relativistic B-W
      // equivalent photon energies
      const double k = 0.5 * ( w2-mp2 )/mp;
      const double kcm = 0.5 * ( w2-mp2 )/w;

      const double epicm = 0.5 * ( w2+mpi2-mp2 )/w, ppicm = sqrt( std::max( 0., epicm*epicm-mpi2 ) );
      const double epi2cm = 0.5 * ( w2 + 4.*mpi2 - mp2 )/w, ppi2cm = sqrt( std::max( 0., epi2cm*epi2cm-4*mpi2 ) );
      const double eetacm = 0.5 * ( w2 + meta2 - mp2 )/w, petacm = sqrt( std::max( 0., eetacm*eetacm-meta2 ) );

      std::array<double,7> width;
      std::array<double,7> pgam;
      std::array<double,7> pwid[3];
      for ( unsigned short i = 0; i < 7; ++i ) {
        const ChristyBostedParameterisation::ResonanceParameters& res = params.resonances[i];

        //--- calculate partial widths
        //----- 1-pion decay mode
        const double x02 = res.x0*res.x0;
        pwid[i][0] = res.width*pow( ppicm/res.ppicmr(), 2.*res.angular_momentum+1. ) * pow( ( res.ppicmr()*res.ppicmr()+x02 )/( ppicm*ppicm+x02 ), res.angular_momentum );
        //----- 2-pion decay mode
        pwid[i][1] = res.width*pow( ppi2cm/res.ppi2cmr(), 2.*res.angular_momentum+4. ) * pow( ( res.ppi2cmr()*res.ppi2cmr()+x02 )/( ppi2cm*ppi2cm+x02 ), res.angular_momentum+2 ) * w / res.mass;
        //----- eta decay mode
        pwid[i][2] = 0.;

        //--- eta decay only for S11's
        if ( i == 1 || i == 4 ) {
          pwid[2][i] = res.width*pow( petacm/res.petacmr(), 2.*res.angular_momentum+1. ) * pow( ( res.petacmr()*res.petacmr()+x02 )/( petacm*petacm+x02 ), res.angular_momentum );
        }

        // virtual photon width
        pgam[i] = res.width * pow( kcm/res.kcmr(), 2 ) * ( res.kcmr()*res.kcmr()+x02 )/( kcm*kcm+x02 );

        width[i] = res.br.singlepi*pwid[i][0] + res.br.doublepi*pwid[i][1] + res.br.eta*pwid[i][2];
      }

      //--- end resonance kinematics and widths calculations

      std::array<double,7> height;
      //--- begin resonance Q^2 dependence calculations
      for ( unsigned short i = 0; i < 7; ++i ) {
        const ChristyBostedParameterisation::ResonanceParameters& res = params.resonances[i];

        if ( sf == 1 ) height[i] = res.A0_T*( 1.+res.fit_parameters[0]*q2/( 1.+res.fit_parameters[1]*q2 ) )/pow( 1.+q2/0.91, res.fit_parameters[2] );
        else           height[i] = res.A0_L*q2/( 1.+res.fit_parameters[3]*q2 )*exp( -q2*res.fit_parameters[4] );
 
        height[i] = height[i]*height[i];
      }

      //--- end resonance Q^2 dependence calculations

      //--- calculate Breit-Wigners for all resonances

      double sig_res = 0.;
      std::array<double,7> sigr;
      for ( unsigned short i = 0; i < 7; ++i ) {
        const ChristyBostedParameterisation::ResonanceParameters res = params.resonances[i];
        const double mass2 = res.mass*res.mass, width2 = res.width*res.width;
        sigr[i] = height[i]*res.kr()/k*res.kcmr()/kcm/res.width * ( res.width*pgam[i] / ( pow( w2 - mass2, 2 ) + mass2*width2 ) );
        //if ( sf == 1 ) sigrsv[i] = sigr[i];
        sig_res += sigr[i];
      }
      sig_res *= w;

      //--- finish resonances / start non-res background calculation
 
      double sig_nr = 0.;
      if ( sf == 1 ) { // transverse
        for ( unsigned short i = 0; i < 2; ++i ) {
          const double wdif = w - ( mp + mpi );
          const double expo = params.continuum.transverse[i].fit_parameters[1]
                            + params.continuum.transverse[i].fit_parameters[2]*q2
                            + params.continuum.transverse[i].fit_parameters[3]*q2*q2;
          sig_nr += params.continuum.transverse[i].sig0 / pow( q2+params.continuum.transverse[i].fit_parameters[0], expo ) * pow( wdif, i+0.5 );
        }

        sig_nr = sig_nr*xpr;
        //sig_nrsv = sig_nr;
      }
      else if ( sf == 2 ) { // longitudinal
        for ( unsigned short i = 0; i < 1; ++i ) {
          const double expo = params.continuum.longitudinal[i].fit_parameters[0]
                            + params.continuum.longitudinal[i].fit_parameters[1];
          sig_nr += params.continuum.longitudinal[i].sig0 * pow( 1.-xpr, expo )/( 1.-xb )
                  * pow( q2/( q2+q20 ), params.continuum.longitudinal[i].fit_parameters[2] )/( q2+q20 )
                  * pow( xpr, params.continuum.longitudinal[i].fit_parameters[3]+params.continuum.longitudinal[i].fit_parameters[4]*t );
        }
      }

      return sig_res + sig_nr;
    }

    ChristyBostedParameterisation
    ChristyBostedParameterisation::standard()
    {
      ChristyBostedParameterisation params;

      params.m0 = 4.2802;
      params.continuum.transverse = {
        ContinuumParameters::DirectionParameters( 246.06, { 0.067469, 1.3501, 0.12054, -0.0038495 } ),
        ContinuumParameters::DirectionParameters( -89.360, { 0.20977, 1.5715, 0.090736, 0.010362 } )
      };
      params.continuum.longitudinal = {
        ContinuumParameters::DirectionParameters( 86.746, { 0., 4.0294, 3.1285, 0.33403, 4.9623 } )
      };

      ResonanceParameters p33; // P33(1232)
      p33.br = ResonanceParameters::BranchingRatios( 1., 0., 1. );
      p33.angular_momentum = 1.;
      //p33.x0 = 0.15;
      p33.x0 = 0.14462;
      p33.mass = 1.2298;
      p33.width = 0.13573;
      p33.fit_parameters = { 4.2291, 1.2598, 2.1242, 19.910, 0.22587 };
      p33.A0_T = 7.7805;
      p33.A0_L = 29.414;
      params.resonances.emplace_back( p33 );

      ResonanceParameters s11_1535; // S11(1535)
      s11_1535.br = ResonanceParameters::BranchingRatios( 0.45, 0.1, 0.45 );
      s11_1535.angular_momentum = 0.;
      s11_1535.x0 = 0.215;
      s11_1535.mass = 1.5304;
      s11_1535.width = 0.220;
      s11_1535.fit_parameters = { 6823.2, 33521., 2.5686, -999., -999. };
      s11_1535.A0_T = 6.3351;
      s11_1535.A0_L = 0.;
      params.resonances.emplace_back( s11_1535 );

      ResonanceParameters d13; // D13(1520)
      d13.br = ResonanceParameters::BranchingRatios( 0.65, 0.35, 0. );
      d13.angular_momentum = 2.;
      d13.x0 = 0.215;
      d13.mass = 1.5057;
      d13.width = 0.082956;
      d13.fit_parameters = { 21.240, 0.055746, 2.4886, 97.046, 0.31042 };
      d13.A0_T = 0.60347;
      d13.A0_L = 157.92;
      params.resonances.emplace_back( d13 );

      ResonanceParameters f15; // F15(1680)
      f15.br = ResonanceParameters::BranchingRatios( 0.65, 0.35, 0. );
      f15.angular_momentum = 3.;
      f15.x0 = 0.215;
      f15.mass = 1.6980;
      f15.width = 0.095782;
      f15.fit_parameters = { -0.28789, 0.18607, 0.063534, 0.038200, 1.2182 };
      f15.A0_T = 2.3305;
      f15.A0_L = 4.2160;
      params.resonances.emplace_back( f15 );

      ResonanceParameters s11_1650; // S11(1650)
      s11_1650.br = ResonanceParameters::BranchingRatios( 0.4, 0.5, 0.1 );
      s11_1650.angular_momentum = 0.;
      s11_1650.x0 = 0.215;
      s11_1650.mass = 1.6650;
      s11_1650.width = 0.10936;
      s11_1650.fit_parameters = { -0.56175, 0.38964, 0.54883, 0.31393, 2.9997 };
      s11_1650.A0_T = 1.9790;
      s11_1650.A0_L = 13.764;
      params.resonances.emplace_back( s11_1650 );

      ResonanceParameters p11; // P11(1440) roper
      p11.br = ResonanceParameters::BranchingRatios( 0.65, 0.35, 0. );
      p11.angular_momentum = 1.;
      p11.x0 = 0.215;
      p11.mass = 1.4333;
      p11.width = 0.37944;
      p11.fit_parameters = { 46.213, 0.19221, 1.9141, 0.053743, 1.3091 };
      p11.A0_T = 0.022506;
      p11.A0_L = 5.5124;
      params.resonances.emplace_back( p11 );

      ResonanceParameters f37; // F37(1950)
      f37.br = ResonanceParameters::BranchingRatios( 0.5, 0.5, 0. );
      f37.angular_momentum = 3.;
      f37.x0 = 0.215;
      f37.mass = 1.9341;
      f37.width = 0.380;
      f37.fit_parameters = { 0., 0., 1., 1.8951, 0.51376 };
      f37.A0_T = 3.4187;
      f37.A0_L = 1.8951;
      params.resonances.emplace_back( f37 );

      return params;
    }

    StructureFunctions
    ChristyBosted( double q2, double xbj )
    {
      const double mp = Particle::massFromPDGId( Particle::Proton ), mp2 = mp*mp;
      const double w2 = mp2 + q2*( 1.-xbj )/xbj;
 
      const double xb = q2/( w2+q2-mp2 );

      const double sigT = resmod507( 1, w2, q2 );
      const double sigL = resmod507( 2, w2, q2 );

      //const double R = sigL/sigT;

      StructureFunctions cb;
      cb.F1 = sigT*(w2-mp2)/8./M_PI/M_PI/Constants::alphaEM/0.3894e3;
      cb.FL = sigL*2.*xb*(w2-mp2)/8./M_PI/M_PI/Constants::alphaEM/0.3894e3;
      return cb;
    }
  }
}

