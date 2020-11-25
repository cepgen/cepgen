//=============================================================================
// NOLI SE TANGERE
#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGen/Core/Exception.h"

#include "CPPProcess.h"

using namespace cepgen;

MadGraphProcess::MadGraphProcess() :
  proc_( new CPPProcess ),
  name_( "XXX_PROC_NAME_XXX" ), descr_( "XXX_PROC_DESCRIPTION_XXX" ),
  central_pdgids_( XXX_OUT_PART_XXX )
{
  CG_INFO( "MadGraphProcess" )
    << "Process considered: " << proc_->name();
}

MadGraphProcess::~MadGraphProcess() = default;

void
MadGraphProcess::initialise( const std::string& param_card )
{
  //--- initialise the process
  try {
    proc_->initProc( param_card );
  } catch ( const char* chr ) {
    throw CG_FATAL( "MadGraphProcess" )
      << "Failed to initialise parameters card at \"" << param_card << "\":\n\t"
      << chr;
  }
  if ( proc_->nprocesses > 1 )
    throw CG_FATAL( "MadGraphProcess" )
      << "Multi-processes matrix elements are not supported!";
  momenta_.clear();
  for ( size_t i = 0; i < proc_->nexternal; ++i )
    momenta_.emplace_back( new double[4]{ 0., 0., 0., proc_->getMasses().at( i ) } );
  if ( proc_->ninitial != 2 )
    throw CG_FATAL( "MadGraphProcess" )
      << "Currently only 2->N processes are supported!";
  particles_.clear();
  particles_.emplace_back( Particle::Role::Parton1, XXX_PART1_XXX, Particle::Status::Incoming );
  particles_.emplace_back( Particle::Role::Parton2, XXX_PART2_XXX, Particle::Status::Incoming );
  for ( size_t i = 0; i < proc_->nexternal; ++i ) {
    particles_.emplace_back( Particle::Role::CentralSystem, central_pdgids_.at( i ), Particle::Status::Undefined );
}

double
MadGraphProcess::eval()
{
  proc_->setMomenta( momenta_ );
  proc_->sigmaKin();
  const double* me = proc_->getMatrixElements();
  if (me[0]!=0.) CG_INFO("")<<me[0];
  //return me[0]*constants::GEVM2_TO_PB;
  return me[0];
}

const Particles&
MadGraphProcess::particles()
{
  const auto& p4 = proc_->getMomenta();
  for ( size_t i = 0; i < p4.size(); ++i )
    particles_[i].setMomentum( p4[i][1], p4[i][2], p4[i][3], p4[i][0] );
  return particles_;
}
//=============================================================================
