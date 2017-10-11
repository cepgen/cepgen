#ifndef CepGen_StructureFunctions_StructureFunctions_h
#define CepGen_StructureFunctions_StructureFunctions_h

namespace CepGen
{
  class StructureFunctions
  {
    public:
      StructureFunctions( double f2=0.0 ) : F2( f2 ), FL( 0.0 ), FM( 0.0 ) {}

      double F2;
      double F1;
      double FL, FM;
  };
}

#endif
