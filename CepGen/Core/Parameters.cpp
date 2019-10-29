#include "CepGen/Parameters.h"

#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/TamingFunction.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/GenericProcess.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/GenericExportHandler.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <iomanip>

namespace cepgen
{
  Parameters::Parameters() :
    general( new ParametersList ),
    store_( false ), total_gen_time_( 0. ), num_gen_events_( 0ul )
  {}

  Parameters::Parameters( Parameters& param ) :
    general( param.general ),
    kinematics( param.kinematics ),
    taming_functions( param.taming_functions ),
    process_( std::move( param.process_ ) ),
    evt_modifiers_( std::move( param.evt_modifiers_ ) ),
    out_module_( std::move( param.out_module_ ) ),
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
    evt_modifiers_ = std::move( param.evt_modifiers_ );
    out_module_ = std::move( param.out_module_ );
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
      << "  structure functions: " << *kinematics.structure_functions;
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

  EventModifier*
  Parameters::eventModifier( size_t i )
  {
    return evt_modifiers_.at( i ).get();
  }

  std::string
  Parameters::eventModifierName( size_t i ) const
  {
    if ( i >= evt_modifiers_.size() )
      return "";
    return evt_modifiers_.at( i )->name();
  }

  void
  Parameters::addModifier( std::unique_ptr<EventModifier> mod )
  {
    evt_modifiers_.emplace_back( std::move( mod ) );
  }

  void
  Parameters::addModifier( EventModifier* mod )
  {
    evt_modifiers_.emplace_back( std::unique_ptr<EventModifier>( mod ) );
  }

  io::GenericExportHandler*
  Parameters::outputModule()
  {
    return out_module_.get();
  }

  void
  Parameters::setOutputModule( std::unique_ptr<io::GenericExportHandler> mod )
  {
    out_module_ = std::move( mod );
  }

  void
  Parameters::setOutputModule( io::GenericExportHandler* mod )
  {
    out_module_.reset( mod );
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
      os << ", " << param->process_->description();
      for ( const auto& par : param->process()->parameters().keys() )
        if ( par != "mode" )
          os << "\n" << std::setw( wt ) << "" << par << ": " << param->process_->parameters().getString( par );
      std::ostringstream proc_mode; proc_mode << param->kinematics.mode;
      if ( param->kinematics.mode != KinematicsMode::invalid )
        os << "\n" << std::setw( wt ) << "Subprocess mode" << ( pretty ? boldify( proc_mode.str().c_str() ) : proc_mode.str() ) << "\n";
    }
    os
      << "\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill( ' ' ) << "\n"
      << std::right << std::setw( wb ) << std::left << std::endl
      << std::setw( wt ) << "Events generation? "
      << ( pretty ? yesno( param->generation_.enabled ) : std::to_string( param->generation_.enabled ) ) << "\n"
      << std::setw( wt ) << "Number of events to generate"
      << ( pretty ? boldify( param->generation_.maxgen ) : std::to_string( param->generation_.maxgen ) ) << "\n";
    if ( param->generation_.num_threads > 1 )
      os
        << std::setw( wt ) << "Number of threads" << param->generation_.num_threads << "\n";
    os
      << std::setw( wt ) << "Number of points to try per bin" << param->generation_.num_points << "\n"
      << std::setw( wt ) << "Integrand treatment"
      << ( pretty ? yesno( param->generation_.treat ) : std::to_string( param->generation_.treat ) ) << "\n"
      << std::setw( wt ) << "Verbosity level " << utils::Logger::get().level << "\n";
    if ( !param->evt_modifiers_.empty() || param->out_module_ || !param->taming_functions.empty() )
      os
        << "\n"
        << std::setfill( '-' ) << std::setw( wb+6 )
        << ( pretty ? boldify( " Event treatment " ) : "Event treatment" ) << std::setfill( ' ' ) << "\n\n";
    if ( !param->evt_modifiers_.empty() ) {
      std::string sep, mod_name = utils::s( "Event modifier", param->evt_modifiers_.size() );
      for ( const auto& mod : param->evt_modifiers_ )
        os
          << std::setw( wt ) << mod_name
          << sep << ( pretty ? boldify( mod->name().c_str() ) : mod->name() ) << "\n", sep = "+ ", mod_name.clear();
      os << "\n";
    }
    if ( param->out_module_ )
      os
        << std::setw( wt ) << "Output module"
        << ( pretty ? boldify( param->out_module_->name().c_str() ) : param->out_module_->name() ) << "\n";
    if ( !param->taming_functions.empty() ) {
      os << std::setw( wt ) << utils::s( "Taming function", param->taming_functions.size() ) << "\n";
      for ( const auto& tf : param->taming_functions )
        os << std::setw( wt ) << ""
          << ( pretty ? boldify( tf.var_orig.c_str() ) : tf.var_orig ) << ": "
          << tf.expr_orig << "\n";
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
    os
      << "\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill( ' ' ) << "\n\n"
      << std::setw( wt ) << "Incoming particles"
      << param->kinematics.incoming_beams.first << ",\n" << std::setw( wt ) << ""
      << param->kinematics.incoming_beams.second << "\n"
      << std::setw( wt ) << "C.m. energy (GeV)" << param->kinematics.sqrtS() << "\n";
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
        os << " * all single " << std::setw( wt-3 ) << PDG::get().name( part_per_lim.first ) << "\n";
        for ( const auto& lim : part_per_lim.second.list() )
          if ( lim.second.valid() )
            os << "   - " << std::setw( wt-5 ) << lim.first << lim.second << "\n";
      }
    }
    os << "\n";
    os << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Proton / remnants " ) : "Proton / remnants" ) << std::setfill( ' ' ) << "\n";
    for ( const auto& lim : param->kinematics.cuts.remnants.list() )
      os << "\n" << std::setw( wt ) << lim.first << lim.second;
    return os
      << "\n"
      << std::setfill('_') << std::setw( wb ) << ""
      << "\n";
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

