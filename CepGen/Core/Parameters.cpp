#include "CepGen/Parameters.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/TamingFunction.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"

namespace CepGen
{
  Parameters::Parameters() :
    general( new ParametersList ),
    hadroniser_max_trials( 5 ),
    taming_functions( new TamingFunctionsCollection ),
    store_( false ), total_gen_time_( 0. ), num_gen_events_( 0 )
  {}

  Parameters::Parameters( Parameters& param ) :
    general( param.general ),
    kinematics( param.kinematics ), integrator( param.integrator ), generation( param.generation ),
    hadroniser_max_trials( param.hadroniser_max_trials ),
    taming_functions( param.taming_functions ),
    process_( std::move( param.process_ ) ),
    hadroniser_( std::move( param.hadroniser_ ) ),
    store_( false ), total_gen_time_( param.total_gen_time_ ), num_gen_events_( param.num_gen_events_ )
  {}

  Parameters::Parameters( const Parameters& param ) :
    general( param.general ),
    kinematics( param.kinematics ), integrator( param.integrator ), generation( param.generation ),
    hadroniser_max_trials( param.hadroniser_max_trials ),
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

  Process::GenericProcess*
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
  Parameters::setProcess( Process::GenericProcess* proc )
  {
    if ( !proc )
      throw CG_FATAL( "Parameters" )
        << "Trying to clone an invalid process!";
    process_.reset( proc );
  }

  void
  Parameters::cloneProcess( const Process::GenericProcess* proc )
  {
    if ( !proc )
      throw CG_FATAL( "Parameters" )
        << "Trying to clone an invalid process!";
    process_ = std::move( proc->clone() );
  }

  Hadroniser::GenericHadroniser*
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
  Parameters::setHadroniser( Hadroniser::GenericHadroniser* hadr )
  {
    hadroniser_.reset( hadr );
  }

  void
  Parameters::dump( std::ostream& out, bool pretty ) const
  {
    std::ostringstream os;

    const int wb = 90, wt = 40;
    os.str( "" );
    os
      << "Parameters dump" << std::left << "\n\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯¯RUN¯INFORMATION¯¯\\_" << std::setfill( ' ' ) << "\n"
      << std::right << std::setw( wb ) << std::left << std::endl
      << std::setw( wt ) << "Process to generate";
    if ( process_ ) {
      os << ( pretty ? boldify( process_->name().c_str() ) : process_->name() ) << "\n"
         << std::setw( wt ) << "" << process_->description();
    }
    else
      os << ( pretty ? boldify( "no process!" ) : "no process!" );
    os
      << "\n"
      << std::setw( wt ) << "Events generation? " << ( pretty ? yesno( generation.enabled ) : std::to_string( generation.enabled ) ) << "\n"
      << std::setw( wt ) << "Number of events to generate" << ( pretty ? boldify( generation.maxgen ) : std::to_string( generation.maxgen ) ) << "\n";
    if ( generation.num_threads > 1 )
      os
        << std::setw( wt ) << "Number of threads" << generation.num_threads << "\n";
    os
      << std::setw( wt ) << "Number of points to try per bin" << generation.num_points << "\n"
      << std::setw( wt ) << "Integrand treatment" << std::boolalpha << generation.treat << "\n"
      << std::setw( wt ) << "Verbosity level " << Logger::get().level << "\n";
    if ( hadroniser_ ) {
      os
        << "\n"
        << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Hadronisation algorithm " ) : "Hadronisation algorithm" ) << std::setfill( ' ' ) << "\n\n"
        << std::setw( wt ) << "Name" << ( pretty ? boldify( hadroniser_->name().c_str() ) : hadroniser_->name() ) << "\n";
    }
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Integration parameters " ) : "Integration parameters" ) << std::setfill( ' ' ) << "\n\n";
    std::ostringstream int_algo; int_algo << integrator.type;
    os
      << std::setw( wt ) << "Integration algorithm" << ( pretty ? boldify( int_algo.str().c_str() ) : int_algo.str() ) << "\n"
      //<< std::setw( wt ) << "Maximum number of iterations" << ( pretty ? boldify( integrator.itvg ) : std::to_string( integrator.itvg ) ) << "\n"
      << std::setw( wt ) << "Number of function calls" << integrator.ncvg << "\n"
      << std::setw( wt ) << "Random number generator seed" << integrator.rng_seed << "\n";
    if ( integrator.rng_engine )
      os
        << std::setw( wt ) << "Random number generator engine" << integrator.rng_engine->name << "\n";
    os
      << "\n"
      << std::setfill('_') << std::setw( wb+3 ) << "_/¯¯EVENTS¯KINEMATICS¯¯\\_" << std::setfill( ' ' ) << "\n\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming particles " ) : "Incoming particles" ) << std::setfill( ' ' ) << "\n\n";
    if ( kinematics.mode != Kinematics::Mode::invalid ) {
      std::ostringstream proc_mode; proc_mode << kinematics.mode;
      os
        << std::setw( wt ) << "Subprocess mode" << ( pretty ? boldify( proc_mode.str().c_str() ) : proc_mode.str() ) << "\n";
      if ( kinematics.mode != Kinematics::Mode::ElasticElastic )
        os << std::setw( wt ) << "Structure functions" << kinematics.structure_functions->type << "\n";
    }
    std::ostringstream ip1, ip2;
    if ( kinematics.incoming_beams.first.hi )
      ip1 << kinematics.incoming_beams.first.hi;
    else
      ip1 << kinematics.incoming_beams.first.pdg;
    if ( kinematics.incoming_beams.second.hi )
      ip2 << kinematics.incoming_beams.second.hi;
    else
      ip2 << kinematics.incoming_beams.second.pdg;
    os
      << std::setw( wt ) << "Incoming particles" << ( pretty ? boldify( ip1.str().c_str() ) : ip1.str() ) << ", " << ( pretty ? boldify( ip2.str().c_str() ) : ip2.str() ) << "\n"
      << std::setw( wt ) << "Momenta (GeV/c)" << kinematics.incoming_beams.first.pz << ", " << kinematics.incoming_beams.second.pz << "\n";
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming partons " ) : "Incoming partons" ) << std::setfill( ' ' ) << "\n\n";
    for ( const auto& lim : kinematics.cuts.initial.list() ) { // map(particles class, limits)
      if ( !lim.second.valid() )
        continue;
      os << std::setw( wt ) << lim.first << lim.second << "\n";
    }
    os
      << "\n"
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing central system " ) : "Outgoing central system" ) << std::setfill( ' ' ) << "\n\n";
    for ( const auto& lim : kinematics.cuts.central.list() ) {
      if ( !lim.second.valid() )
        continue;
      os << std::setw( wt ) << lim.first << lim.second << "\n";
    }
    if ( kinematics.cuts.central_particles.size() > 0 ) {
      os << std::setw( wt ) << ( pretty ? boldify( ">>> per-particle cuts:" ) : ">>> per-particle cuts:" ) << "\n";
      for ( const auto& part_per_lim : kinematics.cuts.central_particles ) {
        os << " * all single " << std::setw( wt-3 ) << part_per_lim.first << "\n";
        for ( const auto& lim : part_per_lim.second.list() ) {
          if ( !lim.second.valid() )
            continue;
          os << "   - " << std::setw( wt-5 ) << lim.first << lim.second << "\n";
        }
      }
    }
    os << "\n";
    os << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Proton / remnants " ) : "Proton / remnants" ) << std::setfill( ' ' ) << "\n\n";
    for ( const auto& lim : kinematics.cuts.remnants.list() )
      os << std::setw( wt ) << lim.first << lim.second << "\n";

    if ( pretty ) {
      CG_INFO( "Parameters" ) << os.str();
    }
    else
      out << os.str();
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
    symmetrise( false ), treat( false ), ngen( 0 ), gen_print_every( 10000 ),
    num_threads( 2 ), num_points( 100 )
  {}
}
