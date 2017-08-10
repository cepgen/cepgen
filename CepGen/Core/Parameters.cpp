#include "CepGen/Parameters.h"

using namespace CepGen;

Parameters::Parameters() :
  process_mode( Kinematics::ElasticElastic ),
  remnant_mode( SuriYennie ),
  in1p( 6500. ), in2p( 6500. ),
  in1pdg( Particle::Proton ), in2pdg( Particle::Proton ),
  pair( Particle::Muon ),
  mcut( Kinematics::BothParticles ),
  minpt( 0. ), maxpt( -1. ),
  minmass( 0. ), maxmass( -1. ),
  minptdiff( 0. ), maxptdiff( -1. ),
  minenergy( 0. ), maxenergy( -1. ),
  mineta( -5. ), maxeta( 5. ),
  minqt( 0. ), maxqt( 500. ),
  minq2( 0. ), maxq2( 1.e5 ),
  minmx( 1.07 ), maxmx( 320. ),
  ncvg( 100000 ), itvg( 10 ), npoints( 100 ), first_run( true ),
  generation( false ), store( false ), maxgen( 0 ),
  symmetrise( true ), ngen( 0 ),
  gpdf( 5 ), spdf( 4 ), qpdf( 12 ),
  hadroniser_max_trials( 5 )
{
  this->last_event = new Event();
  this->file = (std::ofstream*)NULL;
}

Parameters::~Parameters()
{
  delete last_event;
}

void Parameters::setThetaRange( float thetamin, float thetamax )
{
  this->mineta = thetaToEta( thetamax );
  this->maxeta = thetaToEta( thetamin );

  Debugging( Form( "eta(min) = %5.2f => theta(min) = %5.2f"
                   "eta(max) = %5.2f => theta(max) = %5.2f",
                   mineta, thetamin, maxeta, thetamax ) );
}

void Parameters::dump( std::ostream& out, bool pretty ) const
{
  std::ostringstream os;
  os.str( "" ); os << pair; const std::string particles = os.str();
  os.str( "" ); os << mcut; const std::string cutsmode = os.str();

  const int wb = 75, wt = 32;
  os.str( "" );
  os 
    << "Parameters dump" << std::left
    << std::endl << std::endl
    << std::setfill('_') << std::setw( wb ) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill( ' ' ) << std::endl
    << std::right << std::setw( wb ) << std::left << std::endl;
  if ( process )
    os << std::setw( wt ) << "Process to generate" << ( pretty ? boldify( process->name().c_str() ) : process->name() ) << std::endl;
  os
    << std::setw( wt ) << "Events generation? " << ( pretty ? yesno( generation ) : std::to_string( generation ) ) << std::endl
    << std::setw( wt ) << "Number of events to generate" << ( pretty ? boldify( maxgen ) : std::to_string( maxgen ) ) << std::endl
    << std::setw( wt ) << "Events storage? " << yesno( store ) << std::endl
    << std::setw( wt ) << "Verbosity level " << Logger::get().level << std::endl
    << std::setw( wt ) << "Output file opened? " << ( pretty ? yesno( file!=(std::ofstream*)NULL && file->is_open() ) : std::to_string( file!=NULL ) ) << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Vegas integration parameters " ) : "Vegas integration parameters" ) << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setw( wt ) << "Maximum number of iterations" << ( pretty ? boldify( itvg ) : std::to_string( itvg ) ) << std::endl
    << std::setw( wt ) << "Number of function calls" << ncvg << std::endl
    << std::setw( wt ) << "Number of points to try per bin" << npoints << std::endl
    << std::endl
    << std::setfill('_') << std::setw( wb ) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming particles " ) : "Incoming particles" ) << std::setfill( ' ' ) << std::endl
    << std::endl;
  std::ostringstream proc_mode, cut_mode; proc_mode << process_mode; cut_mode << cutsmode;
  std::ostringstream ip1, ip2, op; ip1 << in1pdg; ip2 << in2pdg; op << pair;
  os
    << std::setw( wt ) << "Subprocess mode" << ( pretty ? boldify( proc_mode.str().c_str() ) : proc_mode.str() ) << std::endl
    << std::setw( wt ) << "Incoming particles" << ( pretty ? boldify( ip1.str().c_str() ) : ip1.str() ) << ", " << ( pretty ? boldify( ip2.str().c_str() ) : ip2.str() ) << std::endl
    << std::setw( wt ) << "Momenta (GeV/c)" << in1p << ", " << in2p << std::endl
    << std::setw( wt ) << "Structure functions mode" << remnant_mode << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming partons " ) : "Incoming partons" ) << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setw( wt ) << "Virtuality in range" << ( pretty ? boldify( Form( "%.1f < -t < %.1e", minq2, maxq2 ).c_str() ) : Form( "%.1f < -t < %.1e", minq2, maxq2 ) )<< " GeV**2" << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing leptons " ) : "Outgoing leptons" ) << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setw( wt ) << "Pair" << ( pretty ? boldify( op.str().c_str() ) : op.str() ) << " (" << (int)pair << ")" << std::endl
    << std::setw( wt ) << "Cuts mode" << ( pretty ? boldify( cut_mode.str().c_str() ) : cut_mode.str() ) << std::endl;
  const std::string ptrange = ( maxpt<=0. ) ? Form( ">= %.1f", minpt ) : Form( "%.1f < pT < %.1f", minpt, maxpt ),
                    erange = ( maxenergy<=0. ) ? Form( ">= %.1f", minenergy ) : Form( "%.1f < E < %.1f", minenergy, maxenergy );
  os
    << std::setw( wt ) << "Lepton(s)' pT range" << ( pretty ? boldify( ptrange.c_str() ) : ptrange ) << " GeV/c" << std::endl
    << std::setw( wt ) << "Lepton(s)' energy range" << ( pretty ? boldify( erange.c_str() ) : erange ) << " GeV" << std::endl
    << std::setw( wt ) << "Pseudorapidity range" << ( pretty ? boldify( Form( "%.1f < eta < %.1f", mineta, maxeta ).c_str() ) : Form( "%.1f < eta < %.1f", mineta, maxeta ) ) << std::endl
    //<< std::setw( wt ) << "Polar angle theta in range [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw( 3 ) << maxtheta << "]" << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing remnants" ) : "Outgoing remnants" ) << std::endl
    << std::endl << std::setfill( ' ' );
  if ( hadroniser )
    os << std::setw( wt ) << "Hadronisation algorithm" << ( pretty ? boldify( hadroniser->name().c_str() ) : hadroniser->name() ) << std::endl;
  os << std::setw( wt ) << "Mass range" << ( pretty ? boldify( Form( "%.2f < M(x/y) < %.2f", minmx, maxmx ).c_str() ) : Form( "%.2f < M(x/y) < %.2f", minmx, maxmx ) ) << " GeV/c**2" << std::endl;
  if ( pretty ) { Information( os.str() ); }
  else out << os.str();
}

