#ifndef CepGen_StructureFunctions_SzczurekUleshchenko_h
#define CepGen_StructureFunctions_SzczurekUleshchenko_h

#include "StructureFunctions.h"

extern "C"
{
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& );
  //extern void grv95lo_( double&, double&, double&, double&, double&, double&, double&, double& );
}

namespace CepGen
{
  namespace SF
  {
    class SzczurekUleshchenko : public StructureFunctions
    {
      public:
        SzczurekUleshchenko() {}
        SzczurekUleshchenko operator()( double q2, double xbj ) const;

        double F1;
    };
  }
}

#endif
