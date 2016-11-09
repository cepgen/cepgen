#include "Parameters.h"

Parameters::Parameters() :
  process( 0 ), process_mode( Kinematics::ElasticElastic ),
  remnant_mode( SuriYennie ),
  in1p( 6500. ), in2p( 6500. ),
  in1pdg( Particle::Proton ), in2pdg( Particle::Proton ),
  pair( Particle::Muon ),
  mcut( Kinematics::BothParticles ),
  minpt( 0. ), maxpt( -1. ),
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
  hadroniser( 0 ),
  hadroniser_max_trials( 5 )
{
  this->last_event = new Event();
  this->file = (std::ofstream*)NULL;
}

Parameters::~Parameters()
{
  delete last_event;
}

void Parameters::SetThetaRange( float thetamin_, float thetamax_ )
{
  this->mineta = -log( tan( thetamax_/180.*Constants::Pi/2. ) );
  this->maxeta = -log( tan( thetamin_/180.*Constants::Pi/2. ) );

  Debugging( Form( "eta(min) = %5.2f => theta(min) = %5.2f"
                   "eta(max) = %5.2f => theta(max) = %5.2f",
                   mineta, thetamin_, maxeta, thetamax_ ) );
}

void Parameters::Dump()
{
  std::string cutsmode, particles;
  std::ostringstream os;

  switch (mcut) {
    case 1:           cutsmode = "Vermaseren"; break;
    case 2:           cutsmode = "both leptons"; break;
    case 3:           cutsmode = "single lepton"; break;
    case 0:  default: cutsmode = "none"; break;
  }
  switch (pair) {
    case Particle::Electron:      particles = "electrons"; break;
    case Particle::Muon: default: particles = "muons"; break;
    case Particle::Tau:           particles = "taus"; break;
  }
  const int wb = 75, wt = 32;
  os 
    << "Parameters dump" << std::left
    << std::endl << std::endl
    << std::setfill('_') << std::setw( wb ) << "_/¯ RUN INFORMATION ¯\\_" << std::setfill( ' ' ) << std::endl
    << std::right << std::setw( wb ) << std::left << std::endl;
  if (process)
    os << std::setw( wt ) << "Process to generate" << boldify( process->GetName().c_str() ) << std::endl;
  os
    << std::setw( wt ) << "Events generation? " << yesno( generation ) << std::endl
    << std::setw( wt ) << "Number of events to generate" << boldify( maxgen ) << std::endl
    << std::setw( wt ) << "Events storage? " << yesno( store ) << std::endl
    << std::setw( wt ) << "Debugging mode? " << Logger::GetInstance()->Level << std::endl
    << std::setw( wt ) << "Output file opened? " << yesno( file!=(std::ofstream*)NULL && file->is_open() ) << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb-2 ) << boldify( " Vegas integration parameters " ) << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setw( wt ) << "Maximum number of iterations" << boldify( itvg ) << std::endl
    << std::setw( wt ) << "Number of function calls" << ncvg << std::endl
    << std::setw( wt ) << "Number of points to try per bin" << npoints << std::endl
    << std::endl
    << std::setfill('_') << std::setw( wb ) << "_/¯ EVENTS KINEMATICS ¯\\_" << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb-2 ) << boldify( " Incoming particles " ) << std::setfill( ' ' ) << std::endl
    << std::endl;
  std::ostringstream proc_mode, cut_mode; proc_mode << process_mode; cut_mode << cutsmode;
  std::ostringstream ip1, ip2, op; ip1 << in1pdg; ip2 << in2pdg; op << pair;
  os
    << std::setw( wt ) << "Subprocess mode" << boldify( proc_mode.str().c_str() ) << std::endl
    << std::setw( wt ) << "Incoming particles" << boldify( ip1.str().c_str() ) << ", " << boldify( ip2.str().c_str() ) << std::endl
    << std::setw( wt ) << "Momenta (GeV/c)" << in1p << ", " << in2p << std::endl
    << std::setw( wt ) << "Structure functions mode" << remnant_mode << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb-2 ) << boldify( " Incoming partons " ) << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setw( wt ) << "Virtuality in range" << boldify( Form( "%.1f < -t < %.1e", minq2, maxq2 ).c_str() ) << " GeV**2" << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb-2 ) << boldify( " Outgoing leptons " ) << std::setfill( ' ' ) << std::endl
    << std::endl
    << std::setw( wt ) << "Pair" << boldify( op.str().c_str() ) << " (" << (int)pair << ")" << std::endl
    << std::setw( wt ) << "Cuts mode" << boldify( cut_mode.str().c_str() ) << std::endl;
  const std::string ptrange = ( maxpt<=0. ) ? Form( ">= %.1f", minpt ) : Form( "%.1f < pT < %.1f", minpt, maxpt ),
                    erange = ( maxenergy<=0. ) ? Form( ">= %.1f", minenergy ) : Form( "%.1f < E < %.1f", minenergy, maxenergy );
  os
    << std::setw( wt ) << "Lepton(s)' pT in range" << boldify( ptrange.c_str() ) << " GeV/c" << std::endl
    << std::setw( wt ) << "Lepton(s)' energy in range" << boldify( erange.c_str() ) << " GeV" << std::endl
    << std::setw( wt ) << "Pseudorapidity in range" << boldify( Form( "%.1f < eta < %.1f", mineta, maxeta ).c_str() ) << std::endl
    //<< std::setw( wt ) << "Polar angle theta in range [deg]" << "[" << std::setw(3) << mintheta << ", " << std::setw( 3 ) << maxtheta << "]" << std::endl
    << std::endl
    << std::setfill( '-' ) << std::setw( wb-2 ) << boldify( " Outgoing remnants" ) << std::endl
    << std::endl << std::setfill( ' ' );
  if (hadroniser)
    os << std::setw( wt ) << "Hadronisation algorithm" << boldify( hadroniser->GetName().c_str() ) << std::endl;
  os << std::setw( wt ) << "Mass range" << boldify( Form( "%.2f < M(x/y) < %.2f", minmx, maxmx ).c_str() ) << " GeV/c**2" << std::endl;
  Information( os.str() );
}

bool Parameters::ReadConfigFile(const char* inFile_)
{
  std::ifstream f;
  std::string key, value;
  f.open(inFile_, std::fstream::in);
  if ( !f.is_open() ) {
    return false;
  }

  const unsigned int wdth = 50;
  Debugging( Form( "File '%s' succesfully opened!", inFile_ ) );
  std::ostringstream os;
  os << "Configuration file content :" << "\n";

  os << std::left;
  while (f >> key >> value) {
    //std::cout << std::setw( wdth ) << "[" << key << "] = " << value << std::endl;
    //if ( strncmp( key.c_str(), "#" )==0 ) continue; // FIXME need to ensure there is no extra space before !
    if ( key[0]=='#' ) continue;
    if ( key=="IEND" ) {
      int iend = (int)atoi( value.c_str() );
      if (iend>1) {
        this->generation = true;
      }
    }
    else if ( key=="DEBG" ) {
      Logger::GetInstance()->Level = static_cast<Logger::LoggingLevel>( atoi( value.c_str() ) );
    }
    else if ( key=="NCVG" ) {
      this->ncvg = (int)atoi( value.c_str() );
      os << std::setw( wdth ) << " * Number of function calls:" << this->ncvg << "\n";
    }
    else if ( key=="NCSG" ) {
      this->npoints = (int)atoi( value.c_str() );
      os << std::setw( wdth ) << " * Number of points to probe:" << this->npoints << "\n";
    }
    else if ( key=="ITVG" ) {
      this->itvg = (int)atoi( value.c_str() );
      os << std::setw( wdth ) << " * Number of Vegas iterations:" << this->itvg << "\n";
    }
    else if ( key=="INPP" ) {
      this->in1p = (double)atof( value.c_str() );
      os << std::setw( wdth ) << " * Momentum (1st primary particle):" << this->in1p << " GeV/c\n";
    }
    else if ( key=="INPE" ) {
      this->in2p = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Momentum (2nd primary particle):" << this->in1p << " GeV/c\n";
    }
    else if ( key=="PROC" ) {
      if ( value=="lpair" )       this->process = new GamGamLL;
      else if ( value=="pptoll" ) this->process = new PPtoLL;
      std::ostringstream proc_name; proc_name << this->process;
      os << std::setw( wdth ) << " * Process:" << boldify( proc_name.str() ) << "\n";
    }
    else if ( key=="HADR" ) {
#ifdef PYTHIA6
      if ( value=="pythia6" ) this->hadroniser = new Pythia6Hadroniser;
#endif
#ifdef JETSET
      if ( value=="jetset7" ) this->hadroniser = new Jetset7Hadroniser;
#endif
#ifdef PYTHIA8
      if ( value=="pythia8" ) this->hadroniser = new Pythia8Hadroniser;
#endif
      os << std::setw( wdth ) << " * Hadroniser:" << ( ( this->hadroniser!=0 ) ? this->hadroniser->GetName() : colourise( "*** no hadroniser ***", Colour::Red ) ) << "\n";
    }
    else if ( key=="MODE" ) {
      this->process_mode = static_cast<Kinematics::ProcessMode>( atoi( value.c_str() ) );
      os << std::setw( wdth ) << " * Subprocess' mode:" << static_cast<unsigned int>( this->process_mode ) << " --> " << this->process_mode << "\n";
    }
    else if ( key=="PMOD" or key=="EMOD" ) {
      this->remnant_mode = static_cast<StructureFunctions>( atoi( value.c_str() ) );
      os << std::setw( wdth ) << " * Outgoing primary particles' mode:" << static_cast<unsigned int>( this->remnant_mode )
      	 << " --> " << this->remnant_mode << "\n";
    }
    else if ( key=="PAIR" ) {
      this->pair = static_cast<Particle::ParticleCode>( atoi( value.c_str() ) );
      os << std::setw( wdth ) << " * Outgoing particles' PDG id:" << static_cast<unsigned int>( this->pair )
         << " --> " << this->pair << "\n";
    }
    else if ( key=="MCUT" ) {
      this->mcut = static_cast<Kinematics::Cuts>( atoi( value.c_str() ) );
      os << std::setw( wdth ) << " * Set of cuts to apply on final products:" << this->mcut << "\n";
    }
    else if (key=="PTCT") {
      this->minpt = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Minimal transverse momentum (single central outgoing particle):" << this->minpt << " GeV/c\n";
    }
    else if (key=="ECUT") {
      this->minenergy = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Minimal energy (single central outgoing particle):" << this->minenergy << " GeV\n";
    }
    else if (key=="NGEN") {
      this->maxgen = static_cast<unsigned int>( atoi( value.c_str() ) );
      os << std::setw( wdth ) << " * Number of events to generate:" << boldify( this->maxgen ) << "\n";
    }
    else if (key=="THMN") {
      //this->mintheta = atof( value.c_str() );
      //this->SetThetaRange( atof( value.c_str() ), 0. ); // FIXME FIXME
      os << std::setw( wdth ) << " * Minimal polar production angle for the central particles" << EtaToTheta(mineta) << "\n";
    }
    else if (key=="THMX") {
      //this->maxtheta = atof( value.c_str() );
      //this->SetThetaRange( 0., atof( value.c_str() ) ); //FIXME FIXME
      os << std::setw( wdth ) << " * Maximal polar production angle for the central particles" << EtaToTheta(maxeta) << "\n";
    }
    else if (key=="ETMN") {
      this->mineta = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Minimal pseudo-rapidity (central outgoing particles):" << this->mineta << "\n";
    }
    else if (key=="ETMX") {
      this->maxeta = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Maximal pseudo-rapidity (central outgoing particles):" << this->maxeta << "\n";
    }
    else if (key=="Q2MN") {
      this->minq2 = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Minimal Q^2 (exchanged parton):" << this->minq2 << " GeV^2\n";
    }
    else if (key=="Q2MX") {
      this->maxq2 = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Maximal Q^2 (exchanged parton):" << this->maxq2 << " GeV^2\n";
    }
    else if (key=="MXMN") {
      this->minmx = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Minimal invariant mass of proton remnants:" << this->minmx << " GeV/c^2\n";
    }
    else if (key=="MXMX") {
      this->maxmx = static_cast<float>( atof( value.c_str() ) );
      os << std::setw( wdth ) << " * Maximal invariant mass of proton remnants:" << this->maxmx << " GeV/c^2\n";
    }
    else if (key=="GPDF") {
      this->gpdf = static_cast<unsigned int>( atoi( value.c_str() ) );
      os << std::setw( wdth ) << " * GPDF:" << this->gpdf << "\n";
    }
    else if (key=="SPDF") {
      this->spdf = static_cast<unsigned int>( atoi( value.c_str() ) );
      os << std::setw( wdth ) << " * SPDF:" << this->spdf << "\n";
    }
    else if (key=="QPDF") {
      this->qpdf = static_cast<unsigned int>( atoi( value.c_str() ) );
      os << std::setw( wdth ) << " * QPDF:" << this->qpdf << "\n";
    }
    else {
      InWarning( Form( "Unrecognized argument : [%s] = %s", key.c_str(), value.c_str() ) );
    }
  }
  if ( !this->process ) this->process = new GenericProcess;
  f.close();
  
  Information( os.str() );
  
  return true;
}

bool Parameters::StoreConfigFile( const char* outFile_ )
{
  std::ofstream f;
  f.open( outFile_, std::fstream::out | std::fstream::trunc );
  if ( !f.is_open() ) { return false; }
  // ...
  if ( this->itvg>=0 ) f << "ITVG  " << this->itvg << std::endl;
  if ( this->minenergy!=-1 ) f << "ECUT  " << this->minenergy << std::endl;
  if ( this->minenergy!=-1 ) f << "PTCT  " << this->minpt << std::endl;
  // ...
  f.close();
  return true;
}

