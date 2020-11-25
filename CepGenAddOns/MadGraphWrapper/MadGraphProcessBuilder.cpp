#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

class MadGraphProcessBuilder : public proc::Process
{
  public:
    MadGraphProcessBuilder( const ParametersList& params );
    proc::ProcessPtr clone() const override {
      return proc::ProcessPtr( new MadGraphProcessBuilder( *this ) );
    }
    static std::string description() { return "MadGraph_aMC process builder"; }

    double computeWeight() override;
    void prepareKinematics();
    void fillKinematics( bool ) override;

  private:
    const MadGraphInterface interf_;
    std::shared_ptr<proc::Process> mg5_proc_;
};

MadGraphProcessBuilder::MadGraphProcessBuilder( const ParametersList& params ) :
  Process( params, true ),
  interf_( params )
{
  const auto lib_path = interf_.run();
  loadLibrary( lib_path );
  mg5_proc_ = std::move( proc::ProcessesFactory::get().build( params.get<std::string>( "process" ) ) );
}

void
MadGraphProcessBuilder::prepareKinematics()
{
  mg5_proc_->prepareKinematics();
}

void
MadGraphProcessBuilder::fillKinematics( bool symm )
{
  mg5_proc_->fillKinematics( symm );
}

double
MadGraphProcessBuilder::computeWeight()
{
  return mg5_proc_->computeWeight();
}

REGISTER_PROCESS( "mg5_aMC", MadGraphProcessBuilder )
