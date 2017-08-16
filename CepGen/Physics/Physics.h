#ifndef CepGen_Physics_Physics_h
#define CepGen_Physics_Physics_h

#ifndef __CINT__

#include "Event.h"

extern "C"
{
  //extern void grv95lo_(double&,double&,double&,double&,double&,double&,double&,double&);
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& );
}

/**
 * Lorentz boost of a 4-vector (from CERNLIB)
 * @param pi_ Input 4-vector to boost
 * @param pf_ Output boosted 4-vector
 * @author L. Pape
 * @date 20 Aug 1975
 * @author Ian McLaren (mclareni), CERN/CN
 * @date 14 Feb 1996
 */
void Lorenb( double u_, const CepGen::Particle::Momentum& ps_, double pi_[], double pf_[] );

#endif

#endif
