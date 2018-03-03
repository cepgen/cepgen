#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/Processes/GenericKTProcess.h"

extern "C"
{
  void
  cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
  {
    const CepGen::StructureFunctions::Type mode = ( CepGen::StructureFunctions::Type )sfmode;
    const CepGen::StructureFunctions sf = CepGen::StructureFunctionsBuilder::get( mode, q2, xbj );
    f2 = sf.F2;
    fl = sf.FL;
  }
}

namespace CepGen
{
  namespace Process
  {
    struct FortranKTProcessWrapper : GenericKTProcess
    {
      FortranKTProcessWrapper( const std::string& name, unsigned short parton_pdg, unsigned short outgoing_pdg ) :
        GenericKTProcess( name, std::string( "Fortran wrapped ")+name, 4,
                          { (ParticleCode)parton_pdg, (ParticleCode)parton_pdg },
                          { (ParticleCode)outgoing_pdg, (ParticleCode)outgoing_pdg } ) {}
      ~FortranKTProcessWrapper() {}
      double x[10];
    };
  }
}

extern "C"
{
  extern double test_weight_();
  extern struct {
  } test_params_;
}

#define add_kt_process( name, prefix, parton_pdg, outgoing_pdg ) \
  extern "C" {\
    extern double prefix##_weight_();\
    extern struct {\
    } prefix##_params_;\
  } \

add_kt_process( pptoll_fortran, pptoll, 22, 11 )
