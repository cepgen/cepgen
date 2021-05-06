#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"

#include "CepGen/Processes/Process2to4.h"
#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"

using namespace cepgen;

class MadGraphProcessBuilder : public proc::Process2to4 {
public:
  MadGraphProcessBuilder(const ParametersList& params);
  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new MadGraphProcessBuilder(*this)); }
  static std::string description() { return "MadGraph_aMC process builder"; }

  void prepareProcessKinematics() override;
  double computeCentralMatrixElement() const override;

  /*void preparePhaseSpace() override;
    double computeKTFactorisedMatrixElement() override;
    void fillCentralParticlesKinematics() override;*/

private:
  std::shared_ptr<MadGraphProcess> mg5_proc_;
};

extern std::string madgraph_process_name();

MadGraphProcessBuilder::MadGraphProcessBuilder(const ParametersList& params)
    :  //KTProcess( params, std::array<pdgid_t,2>{}, std::vector<pdgid_t>{} )
      Process2to4(params, std::array<pdgid_t, 2>{}, 0) {
  if (params.has<std::string>("lib"))
    loadLibrary(params.get<std::string>("lib"));
  else {
    const MadGraphInterface interf(params);
    loadLibrary(interf.run());
  }
  //--- once MadGraph process library is loaded into runtime environment
  //    can define its wrapper object
  mg5_proc_.reset(new MadGraphProcess);
}

void
//MadGraphProcessBuilder::preparePhaseSpace()
MadGraphProcessBuilder::prepareProcessKinematics() {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder") << "Process not properly linked!";

  mg5_proc_->initialise(params_.get<std::string>("parametersCard", "param_card.dat"));
}

/*void
MadGraphProcessBuilder::fillCentralParticlesKinematics()
{
  if ( !mg5_proc_ )
    CG_FATAL( "MadGraphProcessBuilder" ) << "Process not properly linked!";

  const auto& parts = mg5_proc_->momenta();
}*/

double
//MadGraphProcessBuilder::computeKTFactorisedMatrixElement()
MadGraphProcessBuilder::computeCentralMatrixElement() const {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder") << "Process not properly linked!";

  //CG_INFO("")<<qt1_*qt1_<<"::"<<t1_<<"\t"<<qt2_*qt2_<<"::"<<t2_;

  // first incoming parton
  /*auto q1 = Momentum::fromPtEtaPhi( qt1_, 0., phi_qt1_ );
  q1.setPz( std::sqrt( t1_-qt1_*qt1_ ) );
  q1.setMass( mg5_proc_->masses()[0] );*/
  //auto q1 = Momentum::fromPtEtaPhiM( qt1_, 0., phi_qt1_, mg5_proc_->masses()[0] );
  mg5_proc_->setMomentum(0, q1_);
  // second incoming parton
  /*auto q2 = Momentum::fromPtEtaPhi( qt2_, 0., phi_qt2_ );
  q2.setPz( std::sqrt( t2_-qt2_*qt2_ ) );
  q2.setMass( mg5_proc_->masses()[1] );*/
  //auto q2 = Momentum::fromPtEtaPhiM( qt2_, 0., phi_qt2_, mg5_proc_->masses()[1] );
  mg5_proc_->setMomentum(1, q2_);

  return mg5_proc_->eval();
}

REGISTER_PROCESS("mg5_aMC", MadGraphProcessBuilder)
