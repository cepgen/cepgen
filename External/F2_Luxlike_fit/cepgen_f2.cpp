#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"

extern "C"
{
  void cepgen_f2_christybosted_( double& xbj, double& q2, double& f2, double& fl )
  {
    CepGen::SF::ChristyBosted cb;
    CepGen::StructureFunctions sf = cb( q2, xbj );
    f2 = sf.F2;
    fl = sf.FL;
  }
  void cepgen_f2_fiorebrasse_( double& xbj, double& q2, double& f2, double& fl )
  {
    CepGen::SF::FioreBrasse fb;
    CepGen::StructureFunctions sf = fb( q2, xbj );
    f2 = sf.F2;
  }
  void cepgen_f2_allm91_( double& xbj, double& q2, double& f2 )
  {
    CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::allm91() );
    CepGen::StructureFunctions sf = allm( q2, xbj );
    f2 = sf.F2;
  }
  void cepgen_f2_allm97_( double& xbj, double& q2, double& f2 )
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

