#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"

extern "C"
{
  void cepgen_f2_christybosted_( double& xbj, double& q2, double& f2, double& fl )
  {
    CepGen::SF::ChristyBosted cb;
    CepGen::StructureFunctions sf = cb( q2, xbj );
    f2 = sf.F2;
    fl = sf.FL;
  }
  void cepgen_f2_allm_( double& xbj, double& q2, double& f2 )
  {
    CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::allm97() );
    CepGen::StructureFunctions sf = allm( q2, xbj );
    f2 = sf.F2;
  }
  void cepgen_f2_gd11p_( double& xbj, double& q2, double& f2 )
  {
    CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::gd11p() );
    CepGen::StructureFunctions sf = allm( q2, xbj );
    f2 = sf.F2;
  }
}

