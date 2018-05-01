#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/Processes/GenericKTProcess.h"

#ifdef __cplusplus
extern "C" {
#endif

  void
  cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
  {
    using namespace CepGen;
    const StructureFunctions::Type mode = ( StructureFunctions::Type )sfmode;
    const StructureFunctions sf = StructureFunctionsBuilder::get( mode, q2, xbj );
    f2 = sf.F2;
    fl = sf.FL;
  }

  double
  cepgen_kt_flux_( int& fmode, double& kt2, double& x, double& mx, int& sfmode )
  {
    using namespace CepGen;
    using namespace CepGen::Process;
    const StructureFunctions::Type sf_mode = ( StructureFunctions::Type )sfmode;
    const GenericKTProcess::Flux f_mode = ( GenericKTProcess::Flux )fmode;
    return GenericKTProcess::flux( f_mode, kt2, x, sf_mode, mx );
  }

  double
  cepgen_kt_flux_hi_( int& fmode, double& kt2, double& x, int& a, int& z )
  {
    using namespace CepGen::Process;
    const GenericKTProcess::Flux f_mode = ( GenericKTProcess::Flux )fmode;
    const GenericKTProcess::HeavyIon hi{ ( unsigned short )a, ( unsigned short )z };
    return GenericKTProcess::flux( f_mode, kt2, x, hi );
  }

#ifdef __cplusplus
}
#endif

