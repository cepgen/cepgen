#include "CepGen/Parameters.h"

#include "CepGen/Integration/Integrator.h"

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/TimeKeeper.h"
#include "CepGen/Utils/String.h"

#include <iomanip>

namespace cepgen
{
  Parameters::Parameters() :
    general( new ParametersList ), integrator( new ParametersList ),
    total_gen_time_( 0. ), num_gen_events_( 0ul ),
    generation_( ParametersList() )
  {}

  Parameters::Parameters( Parameters& param ) :
    general( param.general ), integrator( param.integrator ),
    kinematics( param.kinematics ),
    process_( std::move( param.process_ ) ),
    evt_modifiers_( std::move( param.evt_modifiers_ ) ),
    out_modules_( std::move( param.out_modules_ ) ),
    taming_functions_( std::move( param.taming_functions_ ) ),
    total_gen_time_( param.total_gen_time_ ), num_gen_events_( param.num_gen_events_ ),
    generation_( param.generation_ ),
    tmr_( std::move( param.tmr_ ) )
  {}

  Parameters::Parameters( const Parameters& param ) :
    general( param.general ), integrator( param.integrator ),
    kinematics( param.kinematics ),
    total_gen_time_( param.total_gen_time_ ), num_gen_events_( param.num_gen_events_ ),
    generation_( param.generation_ )
  {}

  Parameters::~Parameters() // required for unique_ptr initialisation!
  {
    CG_DEBUG( "Parameters" ) << "Destructor called.";
  }

  Parameters&
  Parameters::operator=( Parameters param )
  {
    general = param.general;
    integrator = param.integrator;
    kinematics = param.kinematics;
    process_ = std::move( param.process_ );
    evt_modifiers_ = std::move( param.evt_modifiers_ );
    out_modules_ = std::move( param.out_modules_ );
    taming_functions_ = std::move( param.taming_functions_ );
    total_gen_time_ = param.total_gen_time_;
    num_gen_events_ = param.num_gen_events_;
    generation_ = param.generation_;
    tmr_ = std::move( param.tmr_ );
    return *this;
  }

  void
  Parameters::prepareRun()
  {
    if ( tmr_ )
      tmr_->clear();
    CG_TICKER( tmr_.get() );

    //--- first-run preparation
    if ( !process_ || !process_->first_run )
      return;
    CG_DEBUG( "Parameters" ).log( [&]( auto& dbg ) {
      dbg
        << "Run started for " << process_->name() << " process "
        << std::hex << (void*)process_.get() << std::dec << ".\n\t"
        << "Process mode considered: " << kinematics.mode() << "\n\t"
        << "   first beam: " << kinematics.incoming_beams.first << "\n\t"
        << "  second beam: " << kinematics.incoming_beams.second;
      if ( kinematics.structureFunctions() )
        dbg
          << "  structure functions: " << kinematics.structureFunctions();
    } );
    if ( process_->hasEvent() )
      process_->clearEvent();
    //--- clear the run statistics
    total_gen_time_ = 0.;
    num_gen_events_ = 0ul;
    process_->first_run = false;
  }

  void
  Parameters::setTimeKeeper( utils::TimeKeeper* kpr )
  {
    tmr_.reset( kpr );
  }

  void
  Parameters::addGenerationTime( double gen_time )
  {
    total_gen_time_ += gen_time;
    num_gen_events_++;
  }

  proc::Process&
  Parameters::process()
  {
    return *process_.get();
  }

  const proc::Process&
  Parameters::process() const
  {
    return *process_.get();
  }

  std::string
  Parameters::processName() const
  {
    if ( !process_ )
      return "no process";
    return process_->name();
  }

  void
  Parameters::clearProcess()
  {
    process_.release();
  }

  void
  Parameters::setProcess( std::unique_ptr<proc::Process> proc )
  {
    process_ = std::move( proc );
  }

  void
  Parameters::setProcess( proc::Process* proc )
  {
    if ( !proc )
      throw CG_FATAL( "Parameters" )
        << "Trying to clone an invalid process!";
    process_.reset( proc );
  }

  EventModifier&
  Parameters::eventModifier( size_t i )
  {
    return *evt_modifiers_.at( i );
  }

  void
  Parameters::addModifier( std::unique_ptr<EventModifier> mod )
  {
    evt_modifiers_.emplace_back( std::move( mod ) );
    ( *evt_modifiers_.rbegin() )->setParameters( *this );
  }

  void
  Parameters::addModifier( EventModifier* mod )
  {
    evt_modifiers_.emplace_back( std::unique_ptr<EventModifier>( mod ) );
    ( *evt_modifiers_.rbegin() )->setParameters( *this );
  }

  io::ExportModule&
  Parameters::outputModule( size_t i )
  {
    return *out_modules_.at( i );
  }

  void
  Parameters::addOutputModule( std::unique_ptr<io::ExportModule> mod )
  {
    out_modules_.emplace_back( std::move( mod ) );
  }

  void
  Parameters::addOutputModule( io::ExportModule* mod )
  {
    out_modules_.emplace_back( std::unique_ptr<io::ExportModule>( mod ) );
  }

  void
  Parameters::addTamingFunction( std::unique_ptr<utils::Functional> fct )
  {
    taming_functions_.emplace_back( std::move( fct ) );
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
       << ( pretty ? utils::boldify( param->processName() ) : param->processName() );
    if ( param->process_ ) {
      for ( const auto& par : param->process().parameters().keys( false ) )
        if ( par != "mode" )
          os << "\n" << std::setw( wt ) << "" << par << ": " << param->process_->parameters().getString( par );
      std::ostringstream proc_mode; proc_mode << param->kinematics.mode();
      if ( param->kinematics.mode() != mode::Kinematics::invalid )
        os << "\n" << std::setw( wt ) << "Subprocess mode" << ( pretty ? utils::boldify( proc_mode.str() ) : proc_mode.str() ) << "\n";
    }
    os
      << "\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill( ' ' ) << "\n"
      << std::right << std::setw( wb ) << std::left << std::endl
      << std::setw( wt ) << "Events generation? "
      << ( pretty ? utils::yesno( param->generation_.enabled ) : std::to_string( param->generation_.enabled ) ) << "\n"
      << std::setw( wt ) << "Number of events to generate"
      << ( pretty ? utils::boldify( param->generation_.maxgen ) : std::to_string( param->generation_.maxgen ) ) << "\n";
    if ( param->generation_.num_threads > 1 )
      os
        << std::setw( wt ) << "Number of threads" << param->generation_.num_threads << "\n";
    os
      << std::setw( wt ) << "Number of points to try per bin" << param->generation_.num_points << "\n"
      << std::setw( wt ) << "Verbosity level " << utils::Logger::get().level << "\n";
    if ( !param->evt_modifiers_.empty() || param->out_modules_.empty() || !param->taming_functions_.empty() )
      os
        << "\n"
        << std::setfill( '-' ) << std::setw( wb+6 )
        << ( pretty ? utils::boldify( " Event treatment " ) : "Event treatment" ) << std::setfill( ' ' ) << "\n\n";
    if ( !param->evt_modifiers_.empty() ) {
      std::string mod_name = utils::s( "Event modifier", param->evt_modifiers_.size(), false ), sep;
      for ( const auto& mod : param->evt_modifiers_ )
        os
          << std::setw( wt ) << mod_name
          << sep << ( pretty ? utils::boldify( mod->name() ) : mod->name() ) << "\n", sep = "+ ", mod_name.clear();
      os << "\n";
    }
    if ( !param->out_modules_.empty() ) {
      std::string mod_name = utils::s( "Output module", param->out_modules_.size(), false );
      for ( const auto& mod : param->out_modules_ ) {
        os
          << std::setw( wt ) << mod_name
          << ( pretty ? utils::boldify( mod->name() ) : mod->name() ) << "\n", mod_name.clear();
        for ( const auto& par : mod->parameters().keys( false ) )
          os << std::setw( wt ) << "" << par << ": " << mod->parameters().getString( par ) << "\n";
      }
    }
    if ( !param->taming_functions_.empty() ) {
      os << std::setw( wt ) << utils::s( "Taming function", param->taming_functions_.size(), false ) << "\n";
      for ( const auto& tf : param->taming_functions_ )
        os << std::setw( wt ) << ""
          << tf->variables().at( 0 ) << ": " << tf->expression() << "\n";
    }
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 )
      << ( pretty ? utils::boldify( " Integration parameters " ) : "Integration parameters" ) << std::setfill( ' ' ) << "\n\n"
      << std::setw( wt ) << "Integration"
      << ( pretty ? utils::boldify( param->integrator->name<std::string>( "N/A" ) ) : param->integrator->name<std::string>( "N/A" ) ) << "\n";
    for ( const auto& key : param->integrator->keys( false ) )
      os << std::setw( wt ) << "" << key << ": " << param->integrator->getString( key ) << "\n";
    os
      << "\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill( ' ' ) << "\n\n"
      << std::setw( wt ) << "Incoming particles"
      << param->kinematics.incoming_beams.first << ",\n" << std::setw( wt ) << ""
      << param->kinematics.incoming_beams.second << "\n"
      << std::setw( wt ) << "C.m. energy (GeV)" << param->kinematics.sqrtS() << "\n"
      << std::setw( wt ) << "Form factors" << param->kinematics.formFactors() << "\n";
    if ( param->kinematics.mode() != mode::Kinematics::ElasticElastic
      && param->kinematics.structureFunctions() )
      os << std::setw( wt ) << "Structure functions" << param->kinematics.structureFunctions() << "\n";
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? utils::boldify( " Incoming partons " ) : "Incoming partons" ) << std::setfill( ' ' ) << "\n\n";
    const auto& cuts = param->kinematics.cuts;
    for ( const auto& lim : cuts.initial.list() ) // map(particles class, limits)
      if ( lim.limits.valid() )
        os << std::setw( wt ) << lim.description << lim.limits << "\n";
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? utils::boldify( " Outgoing central system " ) : "Outgoing central system" ) << std::setfill( ' ' ) << "\n\n";
    if ( !param->kinematics.minimum_final_state.empty() ) {
      os << std::setw( wt ) << "Minimum final state";
      std::string sep;
      for ( const auto& part : param->kinematics.minimum_final_state )
        os << sep << PDG::get().name( part ), sep = ", ";
      os << "\n";
    }
    for ( const auto& lim : cuts.central.list() )
      if ( lim.limits.valid() )
        os << std::setw( wt ) << lim.description << lim.limits << "\n";
    if ( cuts.central_particles.size() > 0 ) {
      os << std::setw( wt ) << ( pretty ? utils::boldify( ">>> per-particle cuts:" ) : ">>> per-particle cuts:" ) << "\n";
      for ( const auto& part_per_lim : cuts.central_particles ) {
        os << " * all single " << std::setw( wt-3 ) << PDG::get().name( part_per_lim.first ) << "\n";
        for ( const auto& lim : const_cast<cuts::Central&>( part_per_lim.second ).list() )
          if ( lim.limits.valid() )
            os << "   - " << std::setw( wt-5 ) << lim.description << lim.limits << "\n";
      }
    }
    os << "\n";
    os << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? utils::boldify( " Proton / remnants " ) : "Proton / remnants" ) << std::setfill( ' ' ) << "\n";
    for ( const auto& lim : cuts.remnants.list() )
      os << "\n" << std::setw( wt ) << lim.description << lim.limits;
    return os
      << "\n"
      << std::setfill('_') << std::setw( wb ) << ""
      << "\n";
  }

  //-----------------------------------------------------------------------------------------------

  Parameters::Generation::Generation( const ParametersList& params ) :
    enabled( params.get<bool>( "enabled", false ) ),
    maxgen( params.get<int>( "maxgen", 0 ) ),
    symmetrise( params.get<bool>( "symmetrise", false ) ),
    gen_print_every( params.get<int>( "printEvery", 10000 ) ),
    num_threads( params.get<int>( "numThreads", 2 ) ),
    num_points( params.get<int>( "numPoints", 100 ) )
  {}

  Parameters::Generation::Generation( const Generation& rhs ) :
    enabled( rhs.enabled ), maxgen( rhs.maxgen ),
    symmetrise( rhs.symmetrise ), gen_print_every( rhs.gen_print_every ),
    num_threads( rhs.num_threads ), num_points( rhs.num_points )
  {}
}

