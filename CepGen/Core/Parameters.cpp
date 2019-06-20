#include "CepGen/Parameters.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/TamingFunction.h"

#include "CepGen/Physics/PDG.h"

#include "CepGen/Physics/FormFactors.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include <iomanip>

namespace cepgen
{
  Parameters::Parameters() :
    general( new ParametersList ),
    taming_functions( new utils::TamingFunctionsCollection ),
    store_( false ), total_gen_time_( 0. ), num_gen_events_( 0ul )
  {}

  Parameters::Parameters( Parameters& param ) :
    general( param.general ),
    kinematics( param.kinematics ),
    taming_functions( param.taming_functions ),
    process_( std::move( param.process_ ) ),
    hadroniser_( std::move( param.hadroniser_ ) ),
    store_( false ), total_gen_time_( param.total_gen_time_ ), num_gen_events_( param.num_gen_events_ ),
    integration_( param.integration_ ), generation_( param.generation_ )
  {}

  Parameters::Parameters( const Parameters& param ) :
    general( param.general ),
    kinematics( param.kinematics ),
    taming_functions( param.taming_functions ),
    store_( false ), total_gen_time_( param.total_gen_time_ ), num_gen_events_( param.num_gen_events_ ),
    integration_( param.integration_ ), generation_( param.generation_ )
  {}

  Parameters::~Parameters() // required for unique_ptr initialisation!
  {}

  Parameters&
  Parameters::operator=( Parameters param )
  {
    general = param.general;
    kinematics = param.kinematics;
    taming_functions = param.taming_functions;
    process_ = std::move( param.process_ );
    hadroniser_ = std::move( param.hadroniser_ );
    total_gen_time_ = param.total_gen_time_;
    num_gen_events_ = param.num_gen_events_;
    integration_ = param.integration_;
    generation_ = param.generation_;
    return *this;
  }

  void
  Parameters::setThetaRange( float thetamin, float thetamax )
  {
    kinematics.cuts.central.eta_single = {
      Particle::thetaToEta( thetamax ),
      Particle::thetaToEta( thetamin )
    };

    CG_DEBUG( "Parameters" )
      << "eta in range: " << kinematics.cuts.central.eta_single
      << " => theta(min) = " << thetamin << ", theta(max) = " << thetamax << ".";
  }

  void
  Parameters::prepareRun()
  {
    //--- first-run preparation
    if ( !process_ || !process_->first_run )
      return;
    CG_DEBUG( "Parameters" )
      << "Run started for " << process_->name() << " process "
      << "0x" << std::hex << process_.get() << std::dec << ".\n\t"
      << "Process mode considered: " << kinematics.mode << "\n\t"
      << "   first beam: " << kinematics.incoming_beams.first << "\n\t"
      << "  second beam: " << kinematics.incoming_beams.second << "\n\t"
      << "  structure functions: " << kinematics.structure_functions;
    if ( process_->hasEvent() )
      process_->clearEvent();
    //--- clear the run statistics
    total_gen_time_ = 0.;
    num_gen_events_ = 0ul;
    process_->first_run = false;
  }

  void
  Parameters::addGenerationTime( double gen_time )
  {
    total_gen_time_ += gen_time;
    num_gen_events_++;
  }

  proc::GenericProcess*
  Parameters::process()
  {
    return process_.get();
  }

  const proc::GenericProcess*
  Parameters::process() const
  {
    return process_.get();
  }

  std::string
  Parameters::processName() const
  {
    if ( !process_ )
      return "no process";
    return process_->name();
  }

  void
  Parameters::setProcess( std::unique_ptr<proc::GenericProcess> proc )
  {
    process_ = std::move( proc );
  }

  void
  Parameters::setProcess( proc::GenericProcess* proc )
  {
    if ( !proc )
      throw CG_FATAL( "Parameters" )
        << "Trying to clone an invalid process!";
    process_.reset( proc );
  }

  hadr::GenericHadroniser*
  Parameters::hadroniser()
  {
    return hadroniser_.get();
  }

  std::string
  Parameters::hadroniserName() const
  {
    if ( !hadroniser_ )
      return "";
    return hadroniser_->name();
  }

  void
  Parameters::setHadroniser( std::unique_ptr<hadr::GenericHadroniser> hadr )
  {
    hadroniser_ = std::move( hadr );
  }

  void
  Parameters::setHadroniser( hadr::GenericHadroniser* hadr )
  {
    hadroniser_.reset( hadr );
  }

  std::ostream&
  operator<<( std::ostream& os, const Parameters* param )
  {
    const bool pretty = true;

    const int wb = 90, wt = 33;
    os << std::left
       << "\n"
       << std::setfill('_') << std::setw( wb+3 ) << "_/¯ PROCESS INFORMATION ¯\\_" << std::setfill( ' ' ) << "\n"
       << std::right << std::setw( wb ) << std::left << std::endl
       << std::setw( wt ) << "Process to generate"
       << ( pretty ? boldify( param->processName().c_str() ) : param->processName() );
    if ( param->process_ ) {
      os
        << "\n"
        << std::setw( wt ) << "" << param->process_->description();
      for ( const auto& par : param->process()->parameters().keys() )
        os << "\n" << std::setw( wt ) << "" << par << ": " << param->process_->parameters().getString( par ) << "\n";
    }
    os
      << "\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill( ' ' ) << "\n"
      << std::right << std::setw( wb ) << std::left << std::endl
      << std::setw( wt ) << "Events generation? "
      << ( pretty ? yesno( param->generation_.enabled ) : std::to_string( param->generation_.enabled ) ) << "\n"
      << std::setw( wt ) << "Number of events to generate"
      << ( pretty ? boldify( param->generation_.maxgen ) : std::to_string( param->generation_.maxgen ) ) << "\n";
    if ( param->generation().num_threads > 1 )
      os
        << std::setw( wt ) << "Number of threads" << param->generation_.num_threads << "\n";
    os
      << std::setw( wt ) << "Number of points to try per bin" << param->generation_.num_points << "\n"
      << std::setw( wt ) << "Integrand treatment"
      << ( pretty ? yesno( param->generation_.treat ) : std::to_string( param->generation_.treat ) ) << "\n"
      << std::setw( wt ) << "Verbosity level " << utils::Logger::get().level << "\n";
    if ( param->hadroniser_ ) {
      os
        << "\n"
        << std::setfill( '-' ) << std::setw( wb+6 )
        << ( pretty ? boldify( " Hadronisation algorithm " ) : "Hadronisation algorithm" ) << std::setfill( ' ' ) << "\n\n"
        << std::setw( wt ) << "Name"
        << ( pretty ? boldify( param->hadroniser_->name().c_str() ) : param->hadroniser_->name() ) << "\n";
    }
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 )
      << ( pretty ? boldify( " Integration parameters " ) : "Integration parameters" ) << std::setfill( ' ' ) << "\n\n";
    std::ostringstream int_algo; int_algo << param->integration_.type;
    os
      << std::setw( wt ) << "Integration algorithm"
      << ( pretty ? boldify( int_algo.str().c_str() ) : int_algo.str() ) << "\n"
      << std::setw( wt ) << "Number of function calls" << param->integration_.ncvg << "\n"
      << std::setw( wt ) << "Random number generator seed" << param->integration_.rng_seed << "\n";
    if ( param->integration_.rng_engine )
      os
        << std::setw( wt ) << "Random number generator engine"
        << param->integration_.rng_engine->name << "\n";
    std::ostringstream proc_mode; proc_mode << param->kinematics.mode;
    os
      << "\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill( ' ' ) << "\n\n"
      << std::setw( wt ) << "Incoming particles"
      << param->kinematics.incoming_beams.first << ",\n" << std::setw( wt ) << ""
      << param->kinematics.incoming_beams.second << "\n"
      << std::setw( wt ) << "C.m. energy (GeV)" << param->kinematics.sqrtS() << "\n";
    if ( param->kinematics.mode != KinematicsMode::invalid )
      os << std::setw( wt ) << "Subprocess mode" << ( pretty ? boldify( proc_mode.str().c_str() ) : proc_mode.str() ) << "\n";
    if ( param->kinematics.mode != KinematicsMode::ElasticElastic )
      os << std::setw( wt ) << "Structure functions" << *param->kinematics.structure_functions << "\n";
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming partons " ) : "Incoming partons" ) << std::setfill( ' ' ) << "\n\n";
    for ( const auto& lim : param->kinematics.cuts.initial.list() ) // map(particles class, limits)
      if ( lim.second.valid() )
        os << std::setw( wt ) << lim.first << lim.second << "\n";
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing central system " ) : "Outgoing central system" ) << std::setfill( ' ' ) << "\n\n";
    for ( const auto& lim : param->kinematics.cuts.central.list() )
      if ( lim.second.valid() )
        os << std::setw( wt ) << lim.first << lim.second << "\n";
    if ( param->kinematics.cuts.central_particles.size() > 0 ) {
      os << std::setw( wt ) << ( pretty ? boldify( ">>> per-particle cuts:" ) : ">>> per-particle cuts:" ) << "\n";
      for ( const auto& part_per_lim : param->kinematics.cuts.central_particles ) {
        os << " * all single " << std::setw( wt-3 ) << part_per_lim.first << "\n";
        for ( const auto& lim : part_per_lim.second.list() )
          if ( lim.second.valid() )
            os << "   - " << std::setw( wt-5 ) << lim.first << lim.second << "\n";
      }
    }
    os << "\n";
    os << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Proton / remnants " ) : "Proton / remnants" ) << std::setfill( ' ' ) << "\n";
    for ( const auto& lim : param->kinematics.cuts.remnants.list() )
      os << "\n" << std::setw( wt ) << lim.first << lim.second;
    os << "\n";

    std::ostringstream ff_oss, sf_oss, proc_oss;
    std::string sep;
    for ( const auto& mod : ff::FormFactorsHandler::get().modules() )
      ff_oss << sep << mod, sep = ", ";
    sep.clear();
    for ( const auto& mod : strfun::StructureFunctionsHandler::get().modules() )
      sf_oss << sep << (strfun::Type)mod, sep = ", ";
    sep.clear();
    for ( const auto& mod : proc::ProcessesHandler::get().modules() )
      proc_oss << sep << mod, sep = ", ";
    CG_INFO( "Parameters:factories" ) << "Dump of factories"
      << "\n  List of form factors modellings:\n\t" << ff_oss.str()
      << "\n  List of structure functions modellings:\n\t" << sf_oss.str()
      << "\n  List of processes:\n\t" << proc_oss.str();

    return os;
  }

  //-----------------------------------------------------------------------------------------------

  Parameters::Integration::Integration() :
    type( IntegratorType::Vegas ), ncvg( 50000 ),
    rng_seed( 0 ), rng_engine( (gsl_rng_type*)gsl_rng_mt19937 ),
    vegas_chisq_cut( 1.5 ),
    result( -1. ), err_result( -1. )
  {
    const size_t ndof = 10; // random number of dimensions for VEGAS parameters retrieval
    {
      std::shared_ptr<gsl_monte_vegas_state> tmp_state( gsl_monte_vegas_alloc( ndof ), gsl_monte_vegas_free );
      gsl_monte_vegas_params_get( tmp_state.get(), &vegas );
      vegas.iterations = 10;
    } {
      std::shared_ptr<gsl_monte_miser_state> tmp_state( gsl_monte_miser_alloc( ndof ), gsl_monte_miser_free );
      gsl_monte_miser_params_get( tmp_state.get(), &miser );
    }
  }

  Parameters::Integration::Integration( const Integration& rhs ) :
    type( rhs.type ), ncvg( rhs.ncvg ),
    rng_seed( rhs.rng_seed ), rng_engine( rhs.rng_engine ),
    vegas( rhs.vegas ), vegas_chisq_cut( rhs.vegas_chisq_cut ),
    miser( rhs.miser ),
    result( -1. ), err_result( -1. )
  {}

  Parameters::Integration::~Integration()
  {
    //if ( vegas.ostream && vegas.ostream != stdout && vegas.ostream != stderr )
    //  fclose( vegas.ostream );
  }

  //-----------------------------------------------------------------------------------------------

  Parameters::Generation::Generation() :
    enabled( false ), maxgen( 0 ),
    symmetrise( false ), treat( true ), gen_print_every( 10000 ),
    num_threads( 2 ), num_points( 100 )
  {}

  Parameters::Generation::Generation( const Generation& rhs ) :
    enabled( rhs.enabled ), maxgen( rhs.maxgen ),
    symmetrise( rhs.symmetrise ), treat( rhs.treat ), gen_print_every( rhs.gen_print_every ),
    num_threads( rhs.num_threads ), num_points( rhs.num_points )
  {}
}

