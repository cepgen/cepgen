#include "EPA.h"

namespace EPA
{
  double fEPAmax = 0., fYmin = 0., fYmax = 0.;
  double fS = 0., fME2 = 0., fMP2 = 0.;
  double fElDotPr = 0., fEEl = 0.;
  PhotonMode fMode = InvalidMode;
  Particle::Momentum fProton, fElectron;
  PhysicsBoundaries fBoundaries;

  void
  InitialiseEPA( const Particle& el, const Particle& pr, const PhotonMode& mode, const PhysicsBoundaries& b )
  {
    fProton = pr.momentum(), fElectron = el.momentum();
    fMode = mode;
    fBoundaries = b;
  }

  void
  PrepareEPA()
  {
    const double sqs = CMEnergy( fElectron, fProton );
    fS = sqs*sqs;
    fME2 = fElectron.mass2(); fMP2 = fProton.mass2();

    fElDotPr = fElectron.fourProduct( fProton );
    if ( fMode>TransversalLongitudinal ) { fEEl = fElDotPr/fProton.mass(); } // Evaluate photon flux in proton rest frame: set EEL to approx. 50TeV
    else                                 { fEEl = fElectron.energy(); }

    // Calculate Y bounds from [ALI, A. et al. (1987) Heavy quark physics at HERA. Proc. HERA workshop, Hamburg 1987 (ed. R.D. PECCEI), 395-494].
    const double w12 = pow( fBoundaries.wmin,2 )-fProton.mass2();
    // Use trick for quadratic equations for dymin. See [W.H. PRESS et al. (1988): Numerical Recipes in C. Cambridge (Cambridge Univ. Press), p. 156].
    const double ysqr = sqrt( pow( fS-w12, 2 )-4.*w12*fME2 ),
                 dymax_tmp = ( fS+w12+ysqr ) / (2.*( fS+fME2 ) );
    fYmin = std::max( w12/( dymax_tmp*( fS+fME2 ) ), fBoundaries.zmin );
    fYmax = std::min( fBoundaries.zmax, // absolute maximum of y, irrespective of final state
                      std::min( fS/( fS+fME2 ),
                                ( pow( fBoundaries.wmax, 2 )-fMP2+fBoundaries.q2max ) / ( 2*fElDotPr ) ) );

    // Set maximal photon weight for efficient rejection plane
    const double gq2min_init = std::max( pow( fElectron.mass()*fYmin, 2 ) / ( 1.-fYmin ), fBoundaries.q2min ),
                 gq2max_init = std::min( fYmax*fS, fBoundaries.q2max );

    if ( fMode==WeizsackerWilliams ) { fEPAmax = Constants::AlphaReduced*pow( fYmin-2., 2 ); } // WWA approximation
    else {
      // full transversal spectrum (2) or full longitudinal and transversal (3) spectrum
      const double eqe = gq2min_init/fElectron.E2(),
                   emqe2 = pow( fYmin-eqe/4., 2 ),
                   emsqr = ( pow( fYmin*fElDotPr, 2 )+gq2min_init*fMP2 ) / ( fElDotPr*fElDotPr+fME2*fMP2 );
      if ( emsqr<0. ) { InError( Form( "Problem with sqrt(emsqr), %f, at epamax determination", emsqr ) ); return; }

      fEPAmax = Constants::AlphaReduced*fYmin*sqrt( emsqr );
      if ( fMode==Transversal ) { fEPAmax *= ( 2.*( 1.-fYmin )+emqe2+eqe ) / ( emqe2+eqe ); } // Transversal spectrum
      else                      { fEPAmax *= ( 4.*( 1.-fYmin )+emqe2+eqe ) / ( emqe2+eqe ); } // Longitudinal & transversal spectrum
    }
    fEPAmax *= log( fYmax/fYmin )*log( gq2max_init/gq2min_init );

    Debugging( Form( "Y min/max = %f / %f\n\t"
                     "Maximal EPA: %e", fYmin, fYmax, fEPAmax ) );

  }

  bool
  EPA( double x1, double x2, double x3, double* q2, Particle::Momentum* out_ele, Particle::Momentum* out_gam, double* lf )
  {
    // default output
    *q2 = 0.;

    if ( fEPAmax<=0. ) PrepareEPA();

    DebuggingInsideLoop( Form( "EPA max = %f", fEPAmax ) );

    const double y = fYmin*pow( fYmax/fYmin, x1 );
    const double gq2min = std::max( pow( fElectron.mass()*y,2 )/( 1.-y ), fBoundaries.q2min ),
                 gq2max = std::min( y*fS, fBoundaries.q2max ); // calculate actual Q2_min, Q2_max from Y
    // produce Q2 spect. (1/x weighted shape)
    *q2 = gq2min*pow( gq2max/gq2min, x2 );

    /// EPA - WWA spectrum
    double epat, epal;
    if ( fMode==WeizsackerWilliams ) { // WWA approximation
      const double r = Constants::AlphaReduced/( y*( *q2 ) );
      epat = r*( 2.*( 1.-y )*( 1.-fME2*y*y/( ( 1.-y )*( *q2 ) ) )+y*y );
      epal = r*  2.*( 1.-y );
    }
    else {
      const double eqe = ( *q2 )/fEEl,
                   emqe2 = pow( y-eqe/4., 2 ),
                   emsqr = ( pow( y*fElDotPr, 2 )+( *q2 )*fMP2 ) / ( fElDotPr*fElDotPr+fME2*fMP2 );
      if ( emsqr<0. ) { InError( Form( "Problem with sqrt(emsqr), %f, y/Q2 pair rejected", emsqr ) ); return false; }

      const double r = Constants::AlphaReduced/( *q2 )*sqrt( emsqr )/( emqe2+eqe );
      epat = r*( 2.*( 1.-y )+emqe2+eqe );
      if ( fMode==Transversal ) { epal = 0.;            } // Transversal spectrum
      else                      { epal = r*2.*( 1.-y ); } // Longitudinal & transversal spectrum
    }
    double epa = epat+epal;
    *lf = epal/epa; // longitudinal fraction

    // unweight MC
    const double w = sqrt(y*2.*fElDotPr-(*q2)+fMP2),
                 r = ( w>=fBoundaries.wmin and w<=fBoundaries.wmax )
                   ? y*( *q2 )*log( fYmax/fYmin )*log( gq2max/gq2min )
                   : 0.;
    epa *= r; epat *= r; epal *= r;

    if ( epa>fEPAmax ) fEPAmax = epa;

    const double emy = fElectron.energy()*( 1.-y ),
                 exy = fProton.energy()*( *q2 )/fS,
                 eesc = emy+exy;
    const double cthe = ( emy-exy )/eesc;

    if ( fabs(cthe)>1. ) return false;

    const double theta = acos( cthe ), phi = 2.*M_PI*x3;
    *out_ele = Particle::Momentum::fromPThetaPhi( -sqrt( eesc*eesc-fME2 ), theta, phi ); out_ele->setMass( fElectron.mass() );
    *out_gam = fElectron-( *out_ele );
    out_gam->setMass2( -( *q2 ) );

    return true;
  }
}
