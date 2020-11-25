//=============================================================================
// NOLI SE TANGERE
#include "CepGen/Processes/KTProcess.h"
#include "CepGen/Modules/ProcessesFactory.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Core/Exception.h"

#include "CPPProcess.h"

using namespace cepgen;

class MadGraphProcess : public proc::KTProcess
{
  public:
    MadGraphProcess( const ParametersList& params );
    proc::ProcessPtr clone() const override {
      return proc::ProcessPtr( new MadGraphProcess( *this ) );
    }
    static std::string description() { return "XXX_PROC_DESCRIPTION_XXX"; }

  private:
    void preparePhaseSpace() override;
    double computeKTFactorisedMatrixElement() override;
    void fillCentralParticlesKinematics() override;

    CPPProcess proc_;
    const std::string param_card_;
    Momentum qt_1_, qt_2_;
    std::vector<double*> momenta_;
};

MadGraphProcess::MadGraphProcess( const ParametersList& params ) :
  KTProcess( params, std::array<pdgid_t,2>{}, std::vector<pdgid_t>{} ),
  param_card_( params.get<std::string>( "parametersCard", "param_card.dat" ) )
{
  CG_INFO( "MadGraphProcess" )
    << "Process considered: " << proc_.name();
}

void
MadGraphProcess::preparePhaseSpace()
{
  //--- initialise the process
  try {
    proc_.initProc( param_card_ );
  } catch ( const char* chr ) {
    throw CG_FATAL( "MadGraphProcess" )
      << "Failed to initialise parameters card at \"" << param_card_ << "\":\n\t"
      << chr;
  }
  if ( proc_.nprocesses > 1 )
    throw CG_FATAL( "MadGraphProcess" )
      << "Multi-processes matrix elements are not supported!";
  momenta_.clear();
  for ( size_t i = 0; i < proc_.nexternal; ++i )
    momenta_.emplace_back( new double[4]{ 0., 0., 0., proc_.getMasses().at( i ) } );
  setProducedParticles( std::vector<pdgid_t>( proc_.nexternal-proc_.ninitial, PDG::invalid ) );
  //CG_INFO("")<<
  event_->dump();
//  CG_WARNING("")<<momenta_;
//  CG_FATAL("");
}

double
MadGraphProcess::computeKTFactorisedMatrixElement()
{
  // first incoming parton
  qt_1_ = Momentum::fromPtEtaPhi( qt1_, 0., phi_qt1_ );
  qt_1_.setMass( 0. );
  momenta_[0][0] = qt_1_.energy();
  momenta_[0][1] = qt_1_.px();
  momenta_[0][2] = qt_1_.py();
  momenta_[0][3] = qt_1_.pz();
  // second incoming parton
  qt_2_ = Momentum::fromPtEtaPhi( qt2_, 0., phi_qt2_ );
  qt_2_.setMass( 0. );
  momenta_[1][0] = qt_2_.energy();
  momenta_[1][1] = qt_2_.px();
  momenta_[1][2] = qt_2_.py();
  momenta_[1][3] = qt_2_.pz();

  proc_.setMomenta( momenta_ );
  proc_.sigmaKin();
  const double* me = proc_.getMatrixElements();
  CG_INFO("")<<me[0];
  //return me[0]*constants::GEVM2_TO_PB;
  return me[0];
}

void
MadGraphProcess::fillCentralParticlesKinematics()
{
  const auto& p4 = proc_.getMomenta();
  /*CG_WARNING("")<<":::";
  for(const auto& m:mom_filled)
    CG_WARNING("")<<std::vector<double>( m, m+4 );*/
  (*event_)[Particle::Parton1][0].setMomentum( p4[0][1], p4[0][2], p4[0][3], p4[0][0] );
  (*event_)[Particle::Parton2][0].setMomentum( p4[1][1], p4[1][2], p4[1][3], p4[1][0] );
  for ( size_t i = proc_.ninitial; i < p4.size(); ++i )
    (*event_)[Particle::CentralSystem][i-2].setMomentum( p4[i][1], p4[i][2], p4[i][3], p4[i][0] );
  event_->dump();
}

REGISTER_PROCESS( "XXX_PROC_NAME_XXX", MadGraphProcess )
//=============================================================================
