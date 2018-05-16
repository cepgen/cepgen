#include "LpairHandler.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoFF.h"
#include "CepGen/Processes/PPtoWW.h"

#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

#include <fstream>

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for LPAIR input cards

    LpairHandler::LpairHandler( const char* file ) :
      pair_( PDG::invalid )
    {
      std::ifstream f( file, std::fstream::in );
      if ( !f.is_open() )
        throw CG_FATAL( "LpairHandler" ) << "Failed to parse file \"" << file << "%s\".";

      init( &params_ );


      std::unordered_map<std::string, std::string> m_params;
      std::string key, value;
      std::ostringstream os;
      while ( f >> key >> value ) {
        if ( key[0] == '#' ) continue; // FIXME need to ensure there is no extra space before!
        setParameter( key, value );
        m_params.insert( { key, value } );
        if ( getDescription( key ) != "null" )
          os << "\n>> " << key << " = " << std::setw( 15 ) << getParameter( key )
             << " (" << getDescription( key ) << ")";
      }
      f.close();

      if ( proc_name_ == "lpair" )
        params_.setProcess( new Process::GamGamLL );
      else if ( proc_name_ == "pptoll" || proc_name_ == "pptoff" )
        params_.setProcess( new Process::PPtoFF );
      else if ( proc_name_ == "pptoww" )
        params_.setProcess( new Process::PPtoWW );
      else
        throw CG_FATAL( "LpairHandler" ) << "Unrecognised process name: " << proc_name_ << "!";

      if ( integr_type_ == "plain" )
        params_.integrator.type = Integrator::Type::plain;
      else if ( integr_type_ == "Vegas" )
        params_.integrator.type = Integrator::Type::Vegas;
      else if ( integr_type_ == "MISER" )
        params_.integrator.type = Integrator::Type::MISER;
      else if ( integr_type_ != "" )
        throw CG_FATAL( "LpairHandler" ) << "Unrecognized integrator type: " << integr_type_ << "!";

#ifdef PYTHIA8
      if ( hadr_name_ == "pythia8" )
        params_.setHadroniser( new Hadroniser::Pythia8Hadroniser( params_ ) );
#endif

      if ( m_params.count( "IEND" ) )
        setValue<bool>( "IEND", ( std::stoi( m_params["IEND"] ) > 1 ) );

      //--- for LPAIR: specify the lepton pair to be produced
      if ( pair_ != PDG::invalid )
        params_.kinematics.central_system = { pair_, pair_ };

      CG_INFO( "LpairHandler" ) << "File '" << file << "' succesfully opened!\n\t"
        << "The following parameters are set:" << os.str();
    }

    void
    LpairHandler::init( Parameters* params )
    {
      //-------------------------------------------------------------------------------------------
      // Process/integration/hadronisation parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>( "PROC", "Process name to simulate", &proc_name_ );
      registerParameter<std::string>( "ITYP", "Integration algorithm", &integr_type_ );
      registerParameter<std::string>( "HADR", "Hadronisation algorithm", &hadr_name_ );
      registerParameter<std::string>( "KMRG", "KMR grid interpolation path", &params_.kinematics.kmr_grid_path );

      //-------------------------------------------------------------------------------------------
      // General parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<bool>( "IEND", "Generation type", &params->generation.enabled );
      registerParameter<bool>( "NTRT", "Smoothen the integrand", &params->generation.treat );
      registerParameter<unsigned int>( "DEBG", "Debugging verbosity", (unsigned int*)&Logger::get().level );
      registerParameter<unsigned int>( "NCVG", "Number of function calls", &params->integrator.ncvg );
      registerParameter<unsigned int>( "ITVG", "Number of integration iterations", (unsigned int*)&params->integrator.vegas.iterations );
      registerParameter<unsigned int>( "SEED", "Random generator seed", (unsigned int*)&params->integrator.rng_seed );
      registerParameter<unsigned int>( "NTHR", "Number of threads to use for events generation", &params->generation.num_threads );
      registerParameter<unsigned int>( "MODE", "Subprocess' mode", (unsigned int*)&params->kinematics.mode );
      registerParameter<unsigned int>( "PMOD", "Outgoing primary particles' mode", (unsigned int*)&params->kinematics.structure_functions );
      registerParameter<unsigned int>( "EMOD", "Outgoing primary particles' mode", (unsigned int*)&params->kinematics.structure_functions );
      registerParameter<unsigned int>( "PAIR", "Outgoing particles' PDG id", (unsigned int*)&pair_ );
      registerParameter<unsigned int>( "NCSG", "Number of points to probe", &params->generation.num_points );
      registerParameter<unsigned int>( "NGEN", "Number of events to generate", &params->generation.maxgen );
      registerParameter<unsigned int>( "NPRN", "Number of events before printout", &params->generation.gen_print_every );

      //-------------------------------------------------------------------------------------------
      // Process kinematics parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<double>( "INP1", "Momentum (1st primary particle)", &params->kinematics.incoming_beams.first.pz );
      registerParameter<double>( "INP2", "Momentum (2nd primary particle)", &params->kinematics.incoming_beams.second.pz );
      registerParameter<double>( "INPP", "Momentum (1st primary particle)", &params->kinematics.incoming_beams.first.pz );
      registerParameter<double>( "INPE", "Momentum (2nd primary particle)", &params->kinematics.incoming_beams.second.pz );
      registerParameter<double>( "PTCT", "Minimal transverse momentum (single central outgoing particle)", &params->kinematics.cuts.central.pt_single.min() );
      registerParameter<double>( "MSCT", "Minimal central system mass", &params->kinematics.cuts.central.mass_sum.min() );
      registerParameter<double>( "ECUT", "Minimal energy (single central outgoing particle)", &params->kinematics.cuts.central.energy_single.min() );
      registerParameter<double>( "ETMN", "Minimal pseudo-rapidity (central outgoing particles)", &params->kinematics.cuts.central.eta_single.min() );
      registerParameter<double>( "ETMX", "Maximal pseudo-rapidity (central outgoing particles)", &params->kinematics.cuts.central.eta_single.max() );
      registerParameter<double>( "YMIN", "Minimal rapidity (central outgoing particles)", &params->kinematics.cuts.central.rapidity_single.min() );
      registerParameter<double>( "YMAX", "Maximal rapidity (central outgoing particles)", &params->kinematics.cuts.central.rapidity_single.max() );
      registerParameter<double>( "Q2MN", "Minimal Q² = -q² (exchanged parton)", &params->kinematics.cuts.initial.q2.min() );
      registerParameter<double>( "Q2MX", "Maximal Q² = -q² (exchanged parton)", &params->kinematics.cuts.initial.q2.max() );
      registerParameter<double>( "MXMN", "Minimal invariant mass of proton remnants", &params->kinematics.cuts.remnants.mass_single.min() );
      registerParameter<double>( "MXMX", "Maximal invariant mass of proton remnants", &params->kinematics.cuts.remnants.mass_single.max() );
    }

    void
    LpairHandler::store( const char* file )
    {
      std::ofstream f( file, std::fstream::out | std::fstream::trunc );
      if ( !f.is_open() ) {
        CG_ERROR( "LpairHandler" ) << "Failed to open file \"" << file << "%s\" for writing.";
      }
      for ( const auto& it : p_strings_ )
        if ( it.second.value )
          f << it.first << " = " << *it.second.value << "\n";
      for ( const auto& it : p_ints_ )
        if ( it.second.value )
          f << it.first << " = " << *it.second.value << "\n";
      for ( const auto& it : p_doubles_ )
        if ( it.second.value )
          f << it.first << " = " << *it.second.value << "\n";
      for ( const auto& it : p_bools_ )
        if ( it.second.value )
          f << it.first << " = " << *it.second.value << "\n";
      f.close();
    }

    void
    LpairHandler::setParameter( const std::string& key, const std::string& value )
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
      if ( dd != -999. )
        return std::to_string( dd );

      unsigned int ui = getValue<unsigned int>( key.c_str() );
      if ( ui != 999 )
        return std::to_string( ui );

      //if ( out = getValue<bool>( key.c_str() )  );

      return getValue<std::string>( key.c_str() );
    }

    std::string
    LpairHandler::getDescription( std::string key ) const
    {
      if ( p_strings_.count( key ) )
        return p_strings_.find( key )->second.description;
      if ( p_ints_.count( key ) )
        return p_ints_.find( key )->second.description;
      if ( p_doubles_.count( key ) )
        return p_doubles_.find( key )->second.description;
      if ( p_bools_.count( key ) )
        return p_bools_.find( key )->second.description;
      return "null";
    }
  }
}

