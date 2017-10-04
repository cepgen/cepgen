#include "CepGen/StructureFunctions/ALLM.h"

extern "C"
{
  void cepgen_f2_allm_( double& xbj, double& q2, double& f2 )
  {
    CepGen::StructureFunctions sf = CepGen::SF::ALLM( q2, xbj, CepGen::SF::ALLMParameterisation::allm97() );
    f2 = sf.F2;
  }
  void cepgen_f2_gd11p_( double& xbj, double& q2, double& f2 )
  {
    CepGen::StructureFunctions sf = CepGen::SF::ALLM( q2, xbj, CepGen::SF::ALLMParameterisation::gd11p() );
    f2 = sf.F2;
  }
}

