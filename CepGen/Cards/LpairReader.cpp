#include "LpairReader.h"

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for LPAIR input cards

    LpairReader::LpairReader( const char* file )
    {
      std::ifstream f( file, std::fstream::in );
      if ( !f.is_open() ) {
        FatalError( Form( "Failed to parse file \"%s\"", file ) );
        return;
      }

      std::string proc_name, hadr_name;

      registerParameter<std::string>( "PROC", "Process name to simulate", &proc_name );
      registerParameter<std::string>( "HADR", "Hadronisation algorithm to use", &hadr_name );
      registerParameter<bool>( "IEND", "Generation type", &params_.generation.enabled );
      registerParameter<unsigned int>( "DEBG", "Debugging verbosity", (unsigned int*)&Logger::get().level );
      registerParameter<unsigned int>( "NCVG", "Number of function calls", &params_.vegas.ncvg );
      registerParameter<unsigned int>( "NCSG", "Number of points to probe", &params_.vegas.npoints );
      registerParameter<unsigned int>( "ITVG", "Number of Vegas iterations", &params_.vegas.itvg );
      registerParameter<double>( "INPP", "Momentum (1st primary particle)", &params_.kinematics.in1p );
      registerParameter<double>( "INPE", "Momentum (2nd primary particle)", &params_.kinematics.in2p );
      registerParameter<unsigned int>( "MODE", "Subprocess' mode", (unsigned int*)&params_.kinematics.mode );
      registerParameter<unsigned int>( "PMOD", "Outgoing primary particles' mode", (unsigned int*)&params_.remnant_mode );
      registerParameter<unsigned int>( "EMOD", "Outgoing primary particles' mode", (unsigned int*)&params_.remnant_mode );
      registerParameter<unsigned int>( "PAIR", "Outgoing particles' PDG id", (unsigned int*)&params_.kinematics.pair );
      registerParameter<unsigned int>( "MCUT", "Set of cuts to apply on final products", (unsigned int*)&params_.kinematics.cuts_mode );
      registerParameter<double>( "PTCT", "Minimal transverse momentum (single central outgoing particle)", &params_.kinematics.pt_min );
      registerParameter<double>( "MSCT", "Minimal central system mass", &params_.kinematics.mass_min );
      registerParameter<double>( "ECUT", "Minimal energy (single central outgoing particle)", &params_.kinematics.e_min );
      registerParameter<unsigned int>( "NGEN", "Number of events to generate", &params_.generation.maxgen );
      //registerParameter<double>( "THMN", "Minimal polar production angle for the central particles", &params_.kinematics.eta_min );
      //registerParameter<double>( "THMX", "Maximal polar production angle for the central particles", &params_.kinematics.eta_max );
      registerParameter<double>( "ETMN", "Minimal pseudo-rapidity (central outgoing particles)", &params_.kinematics.eta_min );
      registerParameter<double>( "ETMX", "Maximal pseudo-rapidity (central outgoing particles)", &params_.kinematics.eta_max );
      registerParameter<double>( "Q2MN", "Minimal Q^2 (exchanged parton)", &params_.kinematics.q2_min );
      registerParameter<double>( "Q2MX", "Maximal Q^2 (exchanged parton)", &params_.kinematics.q2_max );
      registerParameter<double>( "MXMN", "Minimal invariant mass of proton remnants", &params_.kinematics.mx_min );
      registerParameter<double>( "MXMX", "Maximal invariant mass of proton remnants", &params_.kinematics.mx_max );
      registerParameter<unsigned int>( "GPDF", "GPDF", &params_.pdflib.gpdf );
      registerParameter<unsigned int>( "SPDF", "SPDF", &params_.pdflib.spdf );
      registerParameter<unsigned int>( "QPDF", "QPDF", &params_.pdflib.qpdf );

      Debugging( Form( "File '%s' succesfully opened!", file ) );

      std::map<std::string, std::string> m_params;
      std::string key, value;
      while ( f >> key >> value ) {
        if ( key[0] == '#' ) continue; // FIXME need to ensure there is no extra space before!
        setParameter( key, value );
        m_params.insert( std::pair<std::string,std::string>( key, value ) );
        if ( getDescription( key ) != "null" ) std::cout << "---> " << getDescription( key ) << " = " << getParameter( key ) << std::endl;
      }
      f.close();

      if      ( proc_name == "lpair" )  params_.setProcess( new Process::GamGamLL() );
      else if ( proc_name == "pptoll" ) params_.setProcess( new Process::PPtoLL() );
      else FatalError( Form( "Unrecognised process name: %s", proc_name.c_str() ) );

#ifdef PYTHIA6
      if ( hadr_name == "pythia6" ) params_.setHadroniser( new Hadroniser::Pythia6Hadroniser );
#endif

      if ( m_params.count( "IEND" ) ) setValue<bool>( "IEND", ( std::stoi( m_params["IEND"] ) > 1 ) );
    }

    void
    LpairReader::store( const char* file ) const
    {
      std::ofstream f( file, std::fstream::out | std::fstream::trunc );
      if ( !f.is_open() ) {
        InError( Form( "Failed to open file \"%s\" for writing", file ) );
        return;
      }
      for ( std::vector<Parameter<std::string> >::const_iterator it = p_strings_.begin(); it != p_strings_.end(); ++it ) {
        if ( it->value && !it->value->empty() ) f << it->key << " = " << *it->value << "\n";
      }
      for ( std::vector<Parameter<unsigned int> >::const_iterator it = p_ints_.begin(); it != p_ints_.end(); ++it ) {
        if ( it->value ) f << it->key << " = " << *it->value << "\n";
      }
      for ( std::vector<Parameter<double> >::const_iterator it = p_doubles_.begin(); it != p_doubles_.end(); ++it ) {
        if ( it->value ) f << it->key << " = " << *it->value << "\n";
      }
      f.close();
    }

    void
    LpairReader::setParameter( std::string key, std::string value )
    {
      try { setValue<double>( key.c_str(), std::stod( value ) ); } catch ( std::invalid_argument& ) {}
      try { setValue<unsigned int>( key.c_str(), std::stoi( value ) ); } catch ( std::invalid_argument& ) {}
      //setValue<bool>( key.c_str(), std::stoi( value ) );
      setValue<std::string>( key.c_str(), value );
    }

    std::string
    LpairReader::getParameter( std::string key ) const
    {
      double dd = getValue<double>( key.c_str() );
      if ( dd != -999. ) return std::to_string( dd );

      unsigned int ui = getValue<unsigned int>( key.c_str() );
      if ( ui != 999 ) return std::to_string( ui );

      //if ( out = getValue<bool>( key.c_str() )  );

      return getValue<std::string>( key.c_str() );
    }

    std::string
    LpairReader::getDescription( std::string key ) const
    {
      for ( std::vector<Parameter<std::string> >::const_iterator it = p_strings_.begin(); it != p_strings_.end(); ++it ) { if ( it->key == key ) return it->description; }
      for ( std::vector<Parameter<unsigned int> >::const_iterator it = p_ints_.begin(); it != p_ints_.end(); ++it ) { if ( it->key == key ) return it->description; }
      for ( std::vector<Parameter<double> >::const_iterator it = p_doubles_.begin(); it != p_doubles_.end(); ++it ) { if ( it->key == key ) return it->description; }
      return "null";
    }
  }
}

