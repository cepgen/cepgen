#ifndef __CINT__
#include "Physics.h"

void
Lorenb( double u_, const CepGen::Particle::Momentum& ps_, double pi_[4], double pf_[4] )
{
  double fn;

  if ( ps_.energy()!=u_ ) {
    pf_[3] = ( pi_[3]*ps_.energy() + pi_[2]*ps_.pz() + pi_[1]*ps_.py() + pi_[0]*ps_.px() ) / u_;
    fn = ( pf_[3] + pi_[3] ) / ( ps_.energy() + u_ );
    for ( unsigned int i=0; i<3; i++ ) { pf_[i] = pi_[i]+fn*ps_[i]; }
  }
  else { std::copy( pi_, pi_+4, pf_ ); }
}

#endif
