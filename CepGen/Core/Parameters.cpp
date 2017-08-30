#include "CepGen/Parameters.h"

namespace CepGen
{
  Parameters::Parameters() :
    store_( false )
  {}

  Parameters::Parameters( Parameters& param ) :
    kinematics( param.kinematics ), vegas( param.vegas ), generation( param.generation ),
    taming_functions( param.taming_functions ),
    process_( std::move( param.process_ ) ),
    store_( param.store_ )
  {}

  Parameters::Parameters( const Parameters& param ) :
    kinematics( param.kinematics ), vegas( param.vegas ), generation( param.generation ),
    taming_functions( param.taming_functions ),
    store_( param.store_ )
  {}

  Parameters::~Parameters()
  {}

  void
  Parameters::setThetaRange( float thetamin, float thetamax )
  {
    kinematics.central_cuts[Cuts::eta_single].in( thetaToEta( thetamax ), thetaToEta( thetamin ) );

    if ( Logger::get().level >= Logger::Debug ) {
      std::ostringstream os; os << kinematics.central_cuts[Cuts::eta_single];
      Debugging( Form( "eta in range: %s => theta(min) = %5.2f, theta(max) = %5.2f",
                       os.str().c_str(), thetamin, thetamax ) );
    }
  }

  void
  Parameters::dump( std::ostream& out, bool pretty ) const
  {
    std::ostringstream os;
    os.str( "" ); os << kinematics.cuts_mode; const std::string cutsmode = os.str();

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
      << std::setw( wt ) << "Verbosity level " << Logger::get().level << std::endl
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Vegas integration parameters " ) : "Vegas integration parameters" ) << std::setfill( ' ' ) << std::endl
      << std::endl
      << std::setw( wt ) << "Maximum number of iterations" << ( pretty ? boldify( vegas.itvg ) : std::to_string( vegas.itvg ) ) << std::endl
      << std::setw( wt ) << "Number of function calls" << vegas.ncvg << std::endl
      << std::setw( wt ) << "Number of points to try per bin" << vegas.npoints << std::endl
      << std::endl
      << std::setfill('_') << std::setw( wb ) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill( ' ' ) << std::endl
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming particles " ) : "Incoming particles" ) << std::setfill( ' ' ) << std::endl
      << std::endl;
    std::ostringstream proc_mode, cut_mode; proc_mode << kinematics.mode; cut_mode << cutsmode;
    std::ostringstream ip1, ip2, op; ip1 << kinematics.inpdg.first; ip2 << kinematics.inpdg.second;
    for ( std::vector<Particle::ParticleCode>::const_iterator cp = kinematics.central_system.begin(); cp != kinematics.central_system.end(); ++cp )
      op << ( cp != kinematics.central_system.begin() ? ", " : "" ) << *cp;
    std::ostringstream q2range; q2range << kinematics.initial_cuts.at( Cuts::q2 );
    os
      << std::setw( wt ) << "Subprocess mode" << ( pretty ? boldify( proc_mode.str().c_str() ) : proc_mode.str() ) << std::endl
      << std::setw( wt ) << "Incoming particles" << ( pretty ? boldify( ip1.str().c_str() ) : ip1.str() ) << ", " << ( pretty ? boldify( ip2.str().c_str() ) : ip2.str() ) << std::endl
      << std::setw( wt ) << "Momenta (GeV/c)" << kinematics.inp.first << ", " << kinematics.inp.second << std::endl
      << std::setw( wt ) << "Structure functions used" << kinematics.structure_functions << std::endl
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming partons " ) : "Incoming partons" ) << std::setfill( ' ' ) << std::endl
      << std::endl
      << std::setw( wt ) << "Virtuality range" << ( pretty ? boldify( q2range.str().c_str() ) : q2range.str().c_str() ) << " GeV**2" << std::endl
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing central system " ) : "Outgoing central system" ) << std::setfill( ' ' ) << std::endl
      << std::endl
      << std::setw( wt ) << "Central particles" << ( pretty ? boldify( op.str().c_str() ) : op.str() ) << std::endl
      << std::setw( wt ) << "Cuts mode" << ( pretty ? boldify( cut_mode.str().c_str() ) : cut_mode.str() ) << std::endl;
    for ( std::map<Cuts::Central, Kinematics::Limits>::const_iterator lim = kinematics.central_cuts.begin(); lim != kinematics.central_cuts.end(); ++lim ) {
      os << std::setw( wt ) << lim->first << lim->second << std::endl;
    }
    /*std::ostringstream ptrange; ptrange << kinematics.pt_single_central;
    std::ostringstream erange; erange << kinematics.e_single_central;
    std::ostringstream etarange; etarange << kinematics.eta_single_central;
    std::ostringstream mxrange; mxrange << kinematics.mass_remnants;
    os
      << std::setw( wt ) << "Lepton(s)' pT range" << ( pretty ? boldify( ptrange.str().c_str() ) : ptrange.str().c_str() ) << " GeV/c" << std::endl
      << std::setw( wt ) << "Lepton(s)' energy range" << ( pretty ? boldify( erange.str().c_str() ) : erange.str().c_str() ) << " GeV" << std::endl
      << std::setw( wt ) << "Pseudorapidity range" << ( pretty ? boldify( etarange.str().c_str() )  : etarange.str().c_str() ) << std::endl
      //<< std::setw( wt ) << "Polar angle theta in range [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw( 3 ) << maxtheta << "]" << std::endl
      << std::endl
      << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing remnants" ) : "Outgoing remnants" ) << std::endl
      << std::endl << std::setfill( ' ' );
    os << std::setw( wt ) << "Mass range" << ( pretty ? boldify( mxrange.str().c_str() ) : mxrange.str().c_str() ) << " GeV/c**2" << std::endl;*/
    if ( pretty ) { Information( os.str() ); }
    else out << os.str();
  }
}
