#include "CepGen/Parameters.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Processes/GenericProcess.h"

namespace CepGen
{
  Parameters::Parameters() :
    store_( false )
  {}

  Parameters::Parameters( Parameters& param ) :
    kinematics( param.kinematics ), integrator( param.integrator ), generation( param.generation ),
    hadroniser_max_trials( 1 ),
    taming_functions( param.taming_functions ),
    process_( std::move( param.process_ ) ), hadroniser_( std::move( param.hadroniser_ ) ),
    store_( param.store_ )
  {}

  Parameters::Parameters( const Parameters& param ) :
    kinematics( param.kinematics ), integrator( param.integrator ), generation( param.generation ),
    hadroniser_max_trials( 1 ),
    taming_functions( param.taming_functions ),
    store_( param.store_ )
  {}

  Parameters::~Parameters()
  {}

  void
  Parameters::setThetaRange( float thetamin, float thetamax )
  {
    kinematics.cuts.central[Cuts::eta_single].in( Particle::thetaToEta( thetamax ), Particle::thetaToEta( thetamin ) );

    if ( Logger::get().level >= Logger::Debug ) {
      std::ostringstream os; os << kinematics.cuts.central[Cuts::eta_single];
      Debugging( Form( "eta in range: %s => theta(min) = %5.2f, theta(max) = %5.2f",
                       os.str().c_str(), thetamin, thetamax ) );
    }
  }

  std::string
  Parameters::processName() const
  {
    if ( process_ ) return process_->name();
    return "no process";
  }

  void
  Parameters::dump( std::ostream& out, bool pretty ) const
  {
    std::ostringstream os;

    const int wb = 75, wt = 32;
    os.str( "" );
    os
      << "Parameters dump" << std::left
      << std::endl << std::endl
      << std::setfill('_') << std::setw( wb ) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill( ' ' ) << std::endl
      << std::right << std::setw( wb ) << std::left << std::endl
      << std::setw( wt ) << "Process to generate";
    if ( process_ ) {
      os << ( pretty ? boldify( process_->name().c_str() ) : process_->name() ) << std::endl
         << std::setw( wt ) << "" << process_->description();
    }
    else
      os << ( pretty ? boldify( "no process!" ) : "no process!" );
    os
      << std::endl
      << std::setw( wt ) << "Events generation? " << ( pretty ? yesno( generation.enabled ) : std::to_string( generation.enabled ) ) << std::endl
      << std::setw( wt ) << "Number of events to generate" << ( pretty ? boldify( generation.maxgen ) : std::to_string( generation.maxgen ) ) << std::endl
      << std::setw( wt ) << "Verbosity level " << Logger::get().level << std::endl;
    if ( hadroniser_ ) {
      os
        << std::endl
        << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Hadronisation algorithm " ) : "Hadronisation algorithm" ) << std::setfill( ' ' ) << std::endl
        << std::endl
        << std::setw( wt ) << "Name" << hadroniser_->name() << std::endl;
    }
    os
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Integration parameters " ) : "Integration parameters" ) << std::setfill( ' ' ) << std::endl
      << std::endl;
    std::ostringstream int_algo; int_algo << integrator.type;
    os
      << std::setw( wt ) << "Integration algorithm" << ( pretty ? boldify( int_algo.str().c_str() ) : int_algo.str() ) << std::endl
      << std::setw( wt ) << "Maximum number of iterations" << ( pretty ? boldify( integrator.itvg ) : std::to_string( integrator.itvg ) ) << std::endl
      << std::setw( wt ) << "Number of function calls" << integrator.ncvg << std::endl
      << std::setw( wt ) << "Number of points to try per bin" << integrator.npoints << std::endl
      << std::setw( wt ) << "Random number generator seed" << integrator.seed << std::endl
      << std::endl
      << std::setfill('_') << std::setw( wb ) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill( ' ' ) << std::endl
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming particles " ) : "Incoming particles" ) << std::setfill( ' ' ) << std::endl
      << std::endl;
    std::ostringstream proc_mode; proc_mode << kinematics.mode;
    std::ostringstream ip1, ip2, op; ip1 << kinematics.inpdg.first; ip2 << kinematics.inpdg.second;
    for ( std::vector<ParticleCode>::const_iterator cp = kinematics.central_system.begin(); cp != kinematics.central_system.end(); ++cp )
      op << ( cp != kinematics.central_system.begin() ? ", " : "" ) << *cp;
    std::ostringstream q2range; q2range << kinematics.cuts.initial.at( Cuts::q2 );
    os
      << std::setw( wt ) << "Subprocess mode" << ( pretty ? boldify( proc_mode.str().c_str() ) : proc_mode.str() ) << std::endl
      << std::setw( wt ) << "Incoming particles" << ( pretty ? boldify( ip1.str().c_str() ) : ip1.str() ) << ", " << ( pretty ? boldify( ip2.str().c_str() ) : ip2.str() ) << std::endl
      << std::setw( wt ) << "Momenta (GeV/c)" << kinematics.inp.first << ", " << kinematics.inp.second << std::endl
      << std::setw( wt ) << "Structure functions" << kinematics.structure_functions << std::endl
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming partons " ) : "Incoming partons" ) << std::setfill( ' ' ) << std::endl
      << std::endl
      << std::setw( wt ) << "Virtuality range" << ( pretty ? boldify( q2range.str().c_str() ) : q2range.str().c_str() ) << " GeV**2" << std::endl
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing central system " ) : "Outgoing central system" ) << std::setfill( ' ' ) << std::endl
      << std::endl
      << std::setw( wt ) << "Central particles" << ( pretty ? boldify( op.str().c_str() ) : op.str() ) << std::endl;
    if ( kinematics.cuts.central.size() > 0 ) {
      for ( std::map<Cuts::Central,Kinematics::Limits>::const_iterator lim = kinematics.cuts.central.begin(); lim != kinematics.cuts.central.end(); ++lim ) {
        os << std::setw( wt ) << lim->first << lim->second << std::endl;
      }
    }
    if ( kinematics.cuts.central_particles.size() > 0 ) {
      os << std::setw( wt ) << "> per-particle cuts:" << std::endl;
      for ( std::map<ParticleCode,std::map<Cuts::Central,Kinematics::Limits> >::const_iterator part_lim = kinematics.cuts.central_particles.begin(); part_lim != kinematics.cuts.central_particles.end(); ++part_lim ) {
        os << " * " << std::setw( wt-3 ) << part_lim->first << std::endl;
        for ( std::map<Cuts::Central,Kinematics::Limits>::const_iterator lim = part_lim->second.begin(); lim != part_lim->second.end(); ++lim ) {
          os << "   - " << std::setw( wt-5 ) << lim->first << lim->second << std::endl;
        }
      }
    }
    os << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Proton / remnants " ) : "Proton / remnants" ) << std::setfill( ' ' ) << std::endl;
    for ( std::map<Cuts::Remnants,Kinematics::Limits>::const_iterator lim = kinematics.cuts.remnants.begin(); lim != kinematics.cuts.remnants.end(); ++lim ) {
      os << std::setw( wt ) << lim->first << lim->second << std::endl;
    }
    if ( pretty ) { Information( os.str() ); }
    else out << os.str();
  }
}
