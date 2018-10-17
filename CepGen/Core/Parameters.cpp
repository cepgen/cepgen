#include "CepGen/Parameters.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/TamingFunction.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <iomanip>

namespace cepgen
{
  Parameters::Parameters() :
    general( new ParametersList ),
    taming_functions( new utils::TamingFunctionsCollection ),
    store_( false ), total_gen_time_( 0. ), num_gen_events_( 0 )
  {}

  Parameters::Parameters( Parameters& param ) :
    general( param.general ),
    kinematics( param.kinematics ), integrator( param.integrator ), generation( param.generation ),
    taming_functions( param.taming_functions ),
    process_( std::move( param.process_ ) ),
    hadroniser_( std::move( param.hadroniser_ ) ),
    store_( false ), total_gen_time_( param.total_gen_time_ ), num_gen_events_( param.num_gen_events_ )
  {}

  Parameters::Parameters( const Parameters& param ) :
    general( param.general ),
    kinematics( param.kinematics ), integrator( param.integrator ), generation( param.generation ),
    taming_functions( param.taming_functions ),
    store_( false ), total_gen_time_( param.total_gen_time_ ), num_gen_events_( param.num_gen_events_ )
  {}

  Parameters::~Parameters() // required for unique_ptr initialisation!
  {}

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
  Parameters::clearRunStatistics()
  {
    total_gen_time_ = 0.;
    num_gen_events_ = 0;
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
  operator<<( std::ostream& os, const Parameters* p )
  {
    const bool pretty = true;

    const int wb = 90, wt = 40;
    os
      << "Parameters dump" << std::left << "\n\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯¯RUN¯INFORMATION¯¯\\_" << std::setfill( ' ' ) << "\n"
      << std::right << std::setw( wb ) << std::left << std::endl
      << std::setw( wt ) << "Process to generate";
    if ( p->process_ ) {
      os << ( pretty ? boldify( p->process_->name().c_str() ) : p->process_->name() ) << "\n"
         << std::setw( wt ) << "" << p->process_->description();
    }
    else
      os << ( pretty ? boldify( "no process!" ) : "no process!" );
    os
      << "\n"
      << std::setw( wt ) << "Events generation? "
      << ( pretty ? yesno( p->generation.enabled ) : std::to_string( p->generation.enabled ) ) << "\n"
      << std::setw( wt ) << "Number of events to generate"
      << ( pretty ? boldify( p->generation.maxgen ) : std::to_string( p->generation.maxgen ) ) << "\n";
    if ( p->generation.num_threads > 1 )
      os
        << std::setw( wt ) << "Number of threads" << p->generation.num_threads << "\n";
    os
      << std::setw( wt ) << "Number of points to try per bin" << p->generation.num_points << "\n"
      << std::setw( wt ) << "Integrand treatment"
      << ( pretty ? yesno( p->generation.treat ) : std::to_string( p->generation.treat ) ) << "\n"
      << std::setw( wt ) << "Verbosity level " << utils::Logger::get().level << "\n";
    if ( p->hadroniser_ ) {
      os
        << "\n"
        << std::setfill( '-' ) << std::setw( wb+6 )
        << ( pretty ? boldify( " Hadronisation algorithm " ) : "Hadronisation algorithm" ) << std::setfill( ' ' ) << "\n\n"
        << std::setw( wt ) << "Name"
        << ( pretty ? boldify( p->hadroniser_->name().c_str() ) : p->hadroniser_->name() ) << "\n";
    }
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 )
      << ( pretty ? boldify( " Integration parameters " ) : "Integration parameters" ) << std::setfill( ' ' ) << "\n\n";
    std::ostringstream int_algo; int_algo << p->integrator.type;
    os
      << std::setw( wt ) << "Integration algorithm"
      << ( pretty ? boldify( int_algo.str().c_str() ) : int_algo.str() ) << "\n"
      << std::setw( wt ) << "Number of function calls" << p->integrator.ncvg << "\n"
      << std::setw( wt ) << "Random number generator seed" << p->integrator.rng_seed << "\n";
    if ( p->integrator.rng_engine )
      os
        << std::setw( wt ) << "Random number generator engine"
        << p->integrator.rng_engine->name << "\n";
    std::ostringstream proc_mode; proc_mode << p->kinematics.mode;
    os
      << "\n"
      << std::setfill('_') << std::setw( wb+3 )
      << "_/¯¯EVENTS¯KINEMATICS¯¯\\_" << std::setfill( ' ' ) << "\n\n"
      << std::setw( wt ) << "Incoming particles"
      << p->kinematics.incoming_beams.first << ",\n" << std::setw( wt ) << ""
      << p->kinematics.incoming_beams.second << "\n";
    if ( p->kinematics.mode != KinematicsMode::invalid )
      os << std::setw( wt ) << "Subprocess mode" << ( pretty ? boldify( proc_mode.str().c_str() ) : proc_mode.str() ) << "\n";
    if ( p->kinematics.mode != KinematicsMode::ElasticElastic )
      os << std::setw( wt ) << "Structure functions" << *p->kinematics.structure_functions << "\n";
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming partons " ) : "Incoming partons" ) << std::setfill( ' ' ) << "\n\n";
    for ( const auto& lim : p->kinematics.cuts.initial.list() ) // map(particles class, limits)
      if ( lim.second.valid() )
        os << std::setw( wt ) << lim.first << lim.second << "\n";
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing central system " ) : "Outgoing central system" ) << std::setfill( ' ' ) << "\n\n";
    for ( const auto& lim : p->kinematics.cuts.central.list() )
      if ( lim.second.valid() )
        os << std::setw( wt ) << lim.first << lim.second << "\n";
    if ( p->kinematics.cuts.central_particles.size() > 0 ) {
      os << std::setw( wt ) << ( pretty ? boldify( ">>> per-particle cuts:" ) : ">>> per-particle cuts:" ) << "\n";
      for ( const auto& part_per_lim : p->kinematics.cuts.central_particles ) {
        os << " * all single " << std::setw( wt-3 ) << part_per_lim.first << "\n";
        for ( const auto& lim : part_per_lim.second.list() )
          if ( lim.second.valid() )
            os << "   - " << std::setw( wt-5 ) << lim.first << lim.second << "\n";
      }
    }
    os << "\n";
    os << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Proton / remnants " ) : "Proton / remnants" ) << std::setfill( ' ' ) << "\n\n";
    for ( const auto& lim : p->kinematics.cuts.remnants.list() )
      os << std::setw( wt ) << lim.first << lim.second << "\n";
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const Parameters& p )
  {
    return os << &p;
  }

  //-----------------------------------------------------------------------------------------------

  Parameters::Integration::Integration() :
    type( Integrator::Type::Vegas ), ncvg( 500000 ),
    rng_seed( 0 ), rng_engine( (gsl_rng_type*)gsl_rng_mt19937 ),
    vegas_chisq_cut( 1.5 ),
    result( -1. ), err_result( -1. )
  {
    const size_t ndof = 10;
    {
      std::shared_ptr<gsl_monte_vegas_state> tmp_state( gsl_monte_vegas_alloc( ndof ), gsl_monte_vegas_free );
      gsl_monte_vegas_params_get( tmp_state.get(), &vegas );
    }
    {
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
    symmetrise( false ), treat( true ), ngen( 0 ), gen_print_every( 10000 ),
    num_threads( 2 ), num_points( 100 )
  {}
}
