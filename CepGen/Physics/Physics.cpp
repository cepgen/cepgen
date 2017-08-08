#ifndef __CINT__
#include "Physics.h"

namespace CepGen
{
  PhysicsBoundaries::PhysicsBoundaries() :
    wmin( 20. ), wmax( 0. ),
    q2min( 4. ), q2max( 100. ),
    zmin( 0. ), zmax( 1. )
  {}

  PhysicsBoundaries::~PhysicsBoundaries()
  {}
}

/*double
GetBRFromProcessId(Particle::ParticleCode vmId_)
{
  switch ((VMDecay)vmId_) {
  case RHO_TO_PIPI:         return 1.0;    // rho0->pi+ pi-
  case OMEGA_TO_PIPI:       return 0.0221; // omega->pi+ pi-
  case PHI_TO_KK:           return 0.491;  // phi->K+ K-
  case PHI_TO_KLKS:         return 0.344;  // phi->KL0 KS0 //FIXME FIXME FIXME
  case JPSI_TO_LL:          return 0.0598; // J/psi->l+ l-
  case PSIP_TO_LLX:         return 0.0425; // psi'->l+ l- X
  case UPS1S_TO_LL:         return 0.0250; // Upsilon(1s)->l+ l-
  case UPS2S_TO_LLX:        return 0.0200; // Upsilon(2s)->l+ l- X
  case UPS3S_TO_LLX:        return 0.0217; // Upsilon(3s)->l+ l- X
  case RHO1450_TO_PIPIRHO0: // rho(1450)->pi+ pi- rho0
  case PHI1680_TO_KKBAR: // phi(1680)->K Kbar
///  case RHO_770_0:         return 1.0;    // rho0->pi+ pi-
  case OMEGA_782:       return 0.0221; // omega->pi+ pi-
  case PHI_TO_KK:           return 0.491;  // phi->K+ K-
  case PHI_TO_KLKS:         return 0.344;  // phi->KL0 KS0 //FIXME FIXME FIXME
  case JPSI_TO_LL:          return 0.0598; // J/psi->l+ l-
  case PSIP_TO_LLX:         return 0.0425; // psi'->l+ l- X
  case UPS1S_TO_LL:         return 0.0250; // Upsilon(1s)->l+ l-
  case UPS2S_TO_LLX:        return 0.0200; // Upsilon(2s)->l+ l- X
  case UPS3S_TO_LLX:        return 0.0217; // Upsilon(3s)->l+ l- X
  case RHO1450_TO_PIPIRHO0: // rho(1450)->pi+ pi- rho0
  case PHI1680_TO_KKBAR: // phi(1680)->K Kbar
  default: return -1;
  }
}*/

/*Particles
VMDecayer(Particle part_, GenericHadroniser *had_)
{
  if (part_.status!=1) {
    error.str(""); error << __PRETTY_FUNCTION__ << " ERROR: Particle has status" << part_.status;
    throw std::runtime_error(error.str());
  }
  if (part_.M()<2.5) {
    error.str(""); error << __PRETTY_FUNCTION__ << " ERROR: Particle has a too small mass (" << part_.M() << " GeV)";
    throw std::runtime_error(error.str());
  }
  //if (iret<=-2)...
  if (!had_->Hadronise(&part_)) {

  }
}*/

void
Lorenb( double u_, const CepGen::Particle::Momentum& ps_, double pi_[4], double pf_[4] )
{
  double fn;

  if ( ps_.energy()!=u_ ) {
    pf_[3] = ( pi_[3]*ps_.energy() + pi_[2]*ps_.pz() + pi_[1]*ps_.py() + pi_[0]*ps_.px() ) / u_;
    fn = ( pf_[3] + pi_[3] ) / ( ps_.energy() + u_ );
    for ( unsigned int i=0; i<3; i++ ) { pf_[i] = pi_[i]+fn*ps_.p( i ); }
  }
  else { std::copy( pi_, pi_+4, pf_ ); }
}

#endif
