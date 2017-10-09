#ifndef CepGen_interface_FortranInterface_h
#define CepGen_interface_FortranInterface_h

extern "C"
{
  enum SFmode {
    SuriYennie          = 1,
    SzczurekUleshchenko = 2,
    BlockDurandHa       = 3,
    ALLM91              = 101,
    ALLM97              = 102,
    GD07p               = 103,
    GD11p               = 104,
    FioreBrasse         = 201,
    ChristyBosted       = 202
  };
  void cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl );
}

#endif
