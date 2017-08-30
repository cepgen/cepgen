#include "LpairHandler.h"

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for LPAIR input cards

    LpairHandler::LpairHandler( const char* file )
    {
      std::ifstream f( file, std::fstream::in );
      if ( !f.is_open() ) {
        FatalError( Form( "Failed to parse file \"%s\"", file ) );
        return;
      }
      init( &params_ );

      std::ostringstream os;
      os << Form( "File '%s' succesfully opened! The following parameters are set:\n", file );

      std::map<std::string, std::string> m_params;
      std::string key, value;
      while ( f >> key >> value ) {
        if ( key[0] == '#' ) continue; // FIXME need to ensure there is no extra space before!
        setParameter( key, value );
        m_params.insert( std::pair<std::string,std::string>( key, value ) );
        if ( getDescription( key ) != "null" ) os << ">> " << key << " = " << std::setw( 15 ) << getParameter( key ) << " (" << getDescription( key ) << ")" << std::endl;
      }
      f.close();

      if      ( proc_name_ == "lpair" )  params_.setProcess( new Process::GamGamLL() );
      else if ( proc_name_ == "pptoll" ) params_.setProcess( new Process::PPtoLL() );
      else FatalError( Form( "Unrecognised process name: %s", proc_name_.c_str() ) );

      if ( m_params.count( "IEND" ) ) setValue<bool>( "IEND", ( std::stoi( m_params["IEND"] ) > 1 ) );
      params_.kinematics.central_system = { pair_, pair_ };
      Information( os.str() );
    }

    void
    LpairHandler::init( Parameters* params )
    {
      registerParameter<std::string>( "PROC", "Process name to simulate", &proc_name_ );

      registerParameter<bool>( "IEND", "Generation type", &params->generation.enabled );

      registerParameter<unsigned int>( "DEBG", "Debugging verbosity", (unsigned int*)&Logger::get().level );
      registerParameter<unsigned int>( "NCVG", "Number of function calls", &params->vegas.ncvg );
      registerParameter<unsigned int>( "NCSG", "Number of points to probe", &params->vegas.npoints );
      registerParameter<unsigned int>( "ITVG", "Number of Vegas iterations", &params->vegas.itvg );
      registerParameter<unsigned int>( "MODE", "Subprocess' mode", (unsigned int*)&params->kinematics.mode );
      registerParameter<unsigned int>( "PMOD", "Outgoing primary particles' mode", (unsigned int*)&params->kinematics.structure_functions );
      registerParameter<unsigned int>( "EMOD", "Outgoing primary particles' mode", (unsigned int*)&params->kinematics.structure_functions );
      registerParameter<unsigned int>( "PAIR", "Outgoing particles' PDG id", (unsigned int*)&pair_ );
      registerParameter<unsigned int>( "NGEN", "Number of events to generate", &params->generation.maxgen );

      registerParameter<double>( "INPP", "Momentum (1st primary particle)", &params->kinematics.inp.first );
      registerParameter<double>( "INPE", "Momentum (2nd primary particle)", &params->kinematics.inp.second );
      registerParameter<double>( "PTCT", "Minimal transverse momentum (single central outgoing particle)", &params->kinematics.central_cuts[Cuts::pt_single].min() );
      registerParameter<double>( "MSCT", "Minimal central system mass", &params->kinematics.central_cuts[Cuts::mass_sum].min() );
      registerParameter<double>( "ECUT", "Minimal energy (single central outgoing particle)", &params->kinematics.central_cuts[Cuts::energy_single].min() );
      //registerParameter<double>( "THMN", "Minimal polar production angle for the central particles", &params->kinematics.eta_min );
      //registerParameter<double>( "THMX", "Maximal polar production angle for the central particles", &params->kinematics.eta_max );
      registerParameter<double>( "ETMN", "Minimal pseudo-rapidity (central outgoing particles)", &params->kinematics.central_cuts[Cuts::eta_single].min() );
      registerParameter<double>( "ETMX", "Maximal pseudo-rapidity (central outgoing particles)", &params->kinematics.central_cuts[Cuts::eta_single].max() );
      registerParameter<double>( "Q2MN", "Minimal Q^2 (exchanged parton)", &params->kinematics.initial_cuts[Cuts::q2].min() );
      registerParameter<double>( "Q2MX", "Maximal Q^2 (exchanged parton)", &params->kinematics.initial_cuts[Cuts::q2].max() );
      registerParameter<double>( "MXMN", "Minimal invariant mass of proton remnants", &params->kinematics.remnant_cuts[Cuts::mass].min() );
      registerParameter<double>( "MXMX", "Maximal invariant mass of proton remnants", &params->kinematics.remnant_cuts[Cuts::mass].max() );
    }

    void
    LpairHandler::store( const char* file )
    {
      std::ofstream f( file, std::fstream::out | std::fstream::trunc );
      if ( !f.is_open() ) {
        InError( Form( "Failed to open file \"%s\" for writing", file ) );
        return;
      }
      for ( const auto& it : p_strings_ ) { if ( it.second.value ) f << it.first << " = " << *it.second.value << "\n"; }
      for ( const auto& it : p_ints_ ) { if ( it.second.value ) f << it.first << " = " << *it.second.value << "\n"; }
      for ( const auto& it : p_doubles_ ) { if ( it.second.value ) f << it.first << " = " << *it.second.value << "\n"; }
      for ( const auto& it : p_bools_ ) { if ( it.second.value ) f << it.first << " = " << *it.second.value << "\n"; }
      f.close();
    }

    void
    LpairHandler::setParameter( std::string key, std::string value )
    {
      try { setValue<double>( key.c_str(), std::stod( value ) ); } catch ( std::invalid_argument& ) {}
      try { setValue<unsigned int>( key.c_str(), std::stoi( value ) ); } catch ( std::invalid_argument& ) {}
      //setValue<bool>( key.c_str(), std::stoi( value ) );
      setValue<std::string>( key.c_str(), value );
    }

    std::string
    LpairHandler::getParameter( std::string key ) const
    {
      double dd = getValue<double>( key.c_str() );
      if ( dd != -999. ) return std::to_string( dd );

      unsigned int ui = getValue<unsigned int>( key.c_str() );
      if ( ui != 999 ) return std::to_string( ui );

      //if ( out = getValue<bool>( key.c_str() )  );

      return getValue<std::string>( key.c_str() );
    }

    std::string
    LpairHandler::getDescription( std::string key ) const
    {
      if ( p_strings_.count( key ) ) return p_strings_.find( key )->second.description;
      if ( p_ints_.count( key ) ) return p_ints_.find( key )->second.description;
      if ( p_doubles_.count( key ) ) return p_doubles_.find( key )->second.description;
      if ( p_bools_.count( key ) ) return p_bools_.find( key )->second.description;
      return "null";
    }
  }
}

