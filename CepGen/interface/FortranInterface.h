#ifndef CepGen_interface_FortranInterface_h
#define CepGen_interface_FortranInterface_h

extern "C"
{
  enum SFmode {
    SuriYennie          = 1,
    SzczurekUleshchenko = 2,
    FioreBrasse         = 3,
    ChristyBosted       = 4,
    BlockDurandHa       = 5
  };
  void cepgen_structure_functions_( int sfmode, double q2, double xbj, double& f2, double& fl );
}

#endif
