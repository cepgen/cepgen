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
    std::shared_ptr<MadGraphProcess> mg5_proc_;
};

extern std::string madgraph_process_name();

MadGraphProcessBuilder::MadGraphProcessBuilder( const ParametersList& params ) :
  KTProcess( params, std::array<pdgid_t,2>{}, std::vector<pdgid_t>{} )
{
  if ( params.has<std::string>( "lib" ) )
    loadLibrary( params.get<std::string>( "lib" ) );
  else {
    const MadGraphInterface interf( params );
    loadLibrary( interf.run() );
  }
  //--- once MadGraph process library is loaded into runtime environment
  //    can define its wrapper object
  mg5_proc_.reset( new MadGraphProcess );
}

void
MadGraphProcessBuilder::preparePhaseSpace()
{
  if ( !mg5_proc_ )
    CG_FATAL( "MadGraphProcessBuilder" ) << "Process not properly linked!";
  mg5_proc_->initialise( params_.get<std::string>( "parametersCard", "param_card.dat" ) );
}

void
MadGraphProcessBuilder::fillCentralParticlesKinematics()
{
  if ( !mg5_proc_ )
    CG_FATAL( "MadGraphProcessBuilder" ) << "Process not properly linked!";
  const auto parts = mg5_proc_->particles();
  //mg5_proc_->fillCentralParticlesKinematics();
}

double
MadGraphProcessBuilder::computeKTFactorisedMatrixElement()
{
  if ( !mg5_proc_ )
    CG_FATAL( "MadGraphProcessBuilder" ) << "Process not properly linked!";

  // first incoming parton
  mg5_proc_->setMomentum( 0, Momentum::fromPtEtaPhi( qt1_, 0., phi_qt1_ ) );
  // second incoming parton
  mg5_proc_->setMomentum( 1, Momentum::fromPtEtaPhi( qt2_, 0., phi_qt2_ ) );

  return mg5_proc_->eval();
}

REGISTER_PROCESS( "mg5_aMC", MadGraphProcessBuilder )
