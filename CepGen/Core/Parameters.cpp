#include "CepGen/Parameters.h"

using namespace CepGen;

Parameters::Parameters() :
  process_mode( Kinematics::ElasticElastic ), remnant_mode( SuriYennie ),
  generation( false ), store( false ), maxgen( 0 ),
  last_event( new Event() ),
  symmetrise( true ), ngen( 0 ),
  hadroniser_max_trials( 5 )
{}

Parameters::Parameters( const Parameters& param ) :
  process_mode( param.process_mode ), remnant_mode( param.remnant_mode ),
  kinematics( param.kinematics ), vegas( param.vegas ),
  generation( param.generation ), store( param.store ), maxgen( param.maxgen ),
  last_event( std::move( param.last_event ) ),
  symmetrise( param.symmetrise ), ngen( param.ngen ), pdflib( param.pdflib ),
  hadroniser_max_trials( param.hadroniser_max_trials ),
  process_( std::move( param.process_ ) )
{}

Parameters::~Parameters()
{}

void
Parameters::setThetaRange( float thetamin, float thetamax )
{
  kinematics.eta_min = thetaToEta( thetamax );
  kinematics.eta_max = thetaToEta( thetamin );

  Debugging( Form( "eta(min) = %5.2f => theta(min) = %5.2f"
                   "eta(max) = %5.2f => theta(max) = %5.2f",
                   kinematics.eta_min, thetamin, kinematics.eta_max, thetamax ) );
}

void
Parameters::dump( std::ostream& out, bool pretty ) const
{
  std::ostringstream os;
  os.str( "" ); os << kinematics.pair; const std::string particles = os.str();
  os.str( "" ); os << kinematics.cuts_mode; const std::string cutsmode = os.str();

  const int wb = 75, wt = 32;
  os.str( "" );
  os 
    << "Parameters dump" << std::left
    << std::endl << std::endl
    << std::setfill('_') << std::setw( wb ) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill( ' ' ) << std::endl
    << std::right << std::setw( wb ) << std::left << std::endl
    << std::setw( wt ) << "Process to generate";
  if ( process_ )
    os << ( pretty ? boldify( process_->name().c_str() ) : process_->name() );
  else
    os << ( pretty ? boldify( "no process!" ) : "no process!" );
  os
    << std::endl
    << std::setw( wt ) << "Events generation? " << ( pretty ? yesno( generation ) : std::to_string( generation ) ) << std::endl
    << std::setw( wt ) << "Number of events to generate" << ( pretty ? boldify( maxgen ) : std::to_string( maxgen ) ) << std::endl
    << std::setw( wt ) << "Events storage? " << yesno( store ) << std::endl
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
  std::ostringstream proc_mode, cut_mode; proc_mode << process_mode; cut_mode << cutsmode;
  std::ostringstream ip1, ip2, op; ip1 << kinematics.in1pdg; ip2 << kinematics.in2pdg; op << kinematics.pair;
  os
    << std::setw( wt ) << "Subprocess mode" << ( pretty ? boldify( proc_mode.str().c_str() ) : proc_mode.str() ) << std::endl
    << std::setw( wt ) << "Incoming particles" << ( pretty ? boldify( ip1.str().c_str() ) : ip1.str() ) << ", " << ( pretty ? boldify( ip2.str().c_str() ) : ip2.str() ) << std::endl
    << std::setw( wt ) << "Momenta (GeV/c)" << kinematics.in1p << ", " << kinematics.in2p << std::endl
    << std::setw( wt ) << "Structure functions mode" << remnant_mode << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Incoming partons " ) : "Incoming partons" ) << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setw( wt ) << "Virtuality in range" << ( pretty ? boldify( Form( "%.1f < -t < %.1e", kinematics.q2_min, kinematics.q2_max ).c_str() ) : Form( "%.1f < -t < %.1e", kinematics.q2_min, kinematics.q2_max ) )<< " GeV**2" << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing leptons " ) : "Outgoing leptons" ) << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setw( wt ) << "Pair" << ( pretty ? boldify( op.str().c_str() ) : op.str() ) << " (" << (int)kinematics.pair << ")" << std::endl
    << std::setw( wt ) << "Cuts mode" << ( pretty ? boldify( cut_mode.str().c_str() ) : cut_mode.str() ) << std::endl;
  const std::string ptrange = ( kinematics.pt_max<=0. ) ? Form( ">= %.1f", kinematics.pt_min ) : Form( "%.1f < pT < %.1f", kinematics.pt_min, kinematics.pt_max ),
                    erange = ( kinematics.e_max<=0. ) ? Form( ">= %.1f", kinematics.e_min ) : Form( "%.1f < E < %.1f", kinematics.e_min, kinematics.e_max );
  os
    << std::setw( wt ) << "Lepton(s)' pT range" << ( pretty ? boldify( ptrange.c_str() ) : ptrange ) << " GeV/c" << std::endl
    << std::setw( wt ) << "Lepton(s)' energy range" << ( pretty ? boldify( erange.c_str() ) : erange ) << " GeV" << std::endl
    << std::setw( wt ) << "Pseudorapidity range" << ( pretty ? boldify( Form( "%.1f < eta < %.1f", kinematics.eta_min, kinematics.eta_max ).c_str() ) : Form( "%.1f < eta < %.1f", kinematics.eta_min, kinematics.eta_max ) ) << std::endl
    //<< std::setw( wt ) << "Polar angle theta in range [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw( 3 ) << maxtheta << "]" << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb+6 ) << ( pretty ? boldify( " Outgoing remnants" ) : "Outgoing remnants" ) << std::endl
    << std::endl << std::setfill( ' ' );
  if ( hadroniser_ )
    os << std::setw( wt ) << "Hadronisation algorithm" << ( pretty ? boldify( hadroniser_->name().c_str() ) : hadroniser_->name() ) << std::endl;
  os << std::setw( wt ) << "Mass range" << ( pretty ? boldify( Form( "%.2f < M(x/y) < %.2f", kinematics.mx_min, kinematics.mx_max ).c_str() ) : Form( "%.2f < M(x/y) < %.2f", kinematics.mx_min, kinematics.mx_max ) ) << " GeV/c**2" << std::endl;
  if ( pretty ) { Information( os.str() ); }
  else out << os.str();
}

