#ifndef CepGen_Physics_Physics_h
#define CepGen_Physics_Physics_h

#ifndef __CINT__

#include "Event.h"

extern "C"
{
  //extern void grv95lo_(double&,double&,double&,double&,double&,double&,double&,double&);
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& );
}

#endif

#endif
