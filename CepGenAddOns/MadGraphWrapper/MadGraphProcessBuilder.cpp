#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"

#include "CepGen/Processes/KTProcess.h"
#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

class MadGraphProcessBuilder : public proc::KTProcess
{
  public:
    MadGraphProcessBuilder( const ParametersList& params );
    proc::ProcessPtr clone() const override {
      return proc::ProcessPtr( new MadGraphProcessBuilder( *this ) );
    }
    static std::string description() { return "MadGraph_aMC process builder"; }

    void preparePhaseSpace() override;
    double computeKTFactorisedMatrixElement() override;
    void fillCentralParticlesKinematics() override;

  private:
    std::shared_ptr<proc::KTProcess> mg5_proc_;
};

extern std::string madgraph_process_name();

MadGraphProcessBuilder::MadGraphProcessBuilder( const ParametersList& params ) :
  KTProcess( params, std::array<pdgid_t,2>{}, std::vector<pdgid_t>{} )
{
  auto proc_name = params.get<std::string>( "process" );
  if ( params.has<std::string>( "lib" ) ) {
    loadLibrary( params.get<std::string>( "lib" ) );
    proc_name = madgraph_process_name();
  }
  else {
    const MadGraphInterface interf( params );
    loadLibrary( interf.run() );
  }
  mg5_proc_.reset( dynamic_cast<proc::KTProcess*>( proc::ProcessesFactory::get().build( proc_name ).release() ) );
}

void
MadGraphProcessBuilder::preparePhaseSpace()
{
  if ( !mg5_proc_ )
    CG_FATAL( "MadGraphProcessBuilder" ) << "Process not properly linked!";
  mg5_proc_->preparePhaseSpace();
}

void
MadGraphProcessBuilder::fillCentralParticlesKinematics()
{
  if ( !mg5_proc_ )
    CG_FATAL( "MadGraphProcessBuilder" ) << "Process not properly linked!";
  mg5_proc_->fillCentralParticlesKinematics();
}

double
MadGraphProcessBuilder::computeKTFactorisedMatrixElement()
{
  if ( !mg5_proc_ )
    CG_FATAL( "MadGraphProcessBuilder" ) << "Process not properly linked!";
  return mg5_proc_->computeKTFactorisedMatrixElement();
}

REGISTER_PROCESS( "mg5_aMC", MadGraphProcessBuilder )
