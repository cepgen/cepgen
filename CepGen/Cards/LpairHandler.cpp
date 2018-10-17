#include "CepGen/Cards/LpairHandler.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ProcessesHandler.h"
#include "CepGen/Core/HadronisersHandler.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/GluonGrid.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/StructureFunctions/Partonic.h"

#include <fstream>

namespace cepgen
{
  namespace card
  {
    const int LpairHandler::kInvalid = 99999;

    //----- specialization for LPAIR input cards

    LpairHandler::LpairHandler( const char* file ) :
      proc_params_( new ParametersList ),
      str_fun_( 11 ), hi_1_( { 0, 0 } ), hi_2_( { 0, 0 } )
    {
      std::ifstream f( file, std::fstream::in );
      if ( !f.is_open() )
        throw CG_FATAL( "LpairHandler" ) << "Failed to parse file \"" << file << "%s\".";

      init( &params_ );

      //--- parse all fields
      std::unordered_map<std::string, std::string> m_params;
      std::string key, value;
      std::ostringstream os;
      while ( f >> key >> value ) {
        if ( key[0] == '#' ) // FIXME need to ensure there is no extra space before!
          continue;
        setParameter( key, value );
        m_params.insert( { key, value } );
        if ( getDescription( key ) != "null" )
          os << "\n>> " << key << " = " << std::setw( 15 ) << getParameter( key )
             << " (" << getDescription( key ) << ")";
      }
      f.close();

      //--- parse the process name
      auto proc = cepgen::proc::ProcessesHandler::get().build( proc_name_, *proc_params_ );
      params_.setProcess( std::move( proc ) );

      //--- parse the structure functions code
      const unsigned long kLHAPDFCodeDec = 10000000, kLHAPDFPartDec = 1000000;
      if ( str_fun_ / kLHAPDFCodeDec == 1 ) { // SF from parton
        params_.kinematics.structure_functions = strfun::Parameterisation::build( strfun::Type::Partonic );
        auto sf = dynamic_cast<strfun::Partonic*>( params_.kinematics.structure_functions.get() );
        const unsigned long icode = str_fun_ % kLHAPDFCodeDec;
        sf->params.pdf_code = icode % kLHAPDFPartDec;
        sf->params.mode = (strfun::Partonic::Parameters::Mode)( icode / kLHAPDFPartDec ); // 0, 1, 2
      }
      else
        params_.kinematics.structure_functions = strfun::Parameterisation::build( (strfun::Type)str_fun_ );

      //--- parse the integration algorithm name
      if ( integr_type_ == "plain" )
        params_.integrator.type = Integrator::Type::plain;
      else if ( integr_type_ == "Vegas" )
        params_.integrator.type = Integrator::Type::Vegas;
      else if ( integr_type_ == "MISER" )
        params_.integrator.type = Integrator::Type::MISER;
      else if ( integr_type_ != "" )
        throw CG_FATAL( "LpairHandler" ) << "Unrecognized integrator type: " << integr_type_ << "!";

      //--- parse the hadronisation algorithm name
      auto hadr = cepgen::hadr::HadronisersHandler::get().build( hadr_name_, ParametersList() );
      hadr->setParameters( params_ );
      params_.setHadroniser( std::move( hadr ) );

      if ( m_params.count( "IEND" ) )
        setValue<bool>( "IEND", ( std::stoi( m_params["IEND"] ) > 1 ) );

      if ( m_params.count( "KMRG" ) && !kmr_grid_path_.empty() )
        kmr::GluonGrid::get( kmr_grid_path_.c_str() );

      //--- check if we are dealing with heavy ions for incoming states
      HeavyIon hi1{ hi_1_.first, (Element)hi_1_.second }, hi2{ hi_2_.first, (Element)hi_2_.second };
      if ( hi1 )
        params_.kinematics.incoming_beams.first.pdg = hi1;
      if ( hi2 )
        params_.kinematics.incoming_beams.second.pdg = hi2;

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
      registerParameter<std::string>( "KMRG", "KMR grid interpolation path", &kmr_grid_path_ );

      //-------------------------------------------------------------------------------------------
      // General parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<bool>( "IEND", "Generation type", &params->generation.enabled );
      registerParameter<bool>( "NTRT", "Smoothen the integrand", &params->generation.treat );
      registerParameter<int>( "DEBG", "Debugging verbosity", (int*)&utils::Logger::get().level );
      registerParameter<int>( "NCVG", "Number of function calls", (int*)&params->integrator.ncvg );
      registerParameter<int>( "ITVG", "Number of integration iterations", (int*)&params->integrator.vegas.iterations );
      registerParameter<int>( "SEED", "Random generator seed", (int*)&params->integrator.rng_seed );
      registerParameter<int>( "NTHR", "Number of threads to use for events generation", (int*)&params->generation.num_threads );
      registerParameter<int>( "MODE", "Subprocess' mode", (int*)&params->kinematics.mode );
      registerParameter<int>( "NCSG", "Number of points to probe", (int*)&params->generation.num_points );
      registerParameter<int>( "NGEN", "Number of events to generate", (int*)&params->generation.maxgen );
      registerParameter<int>( "NPRN", "Number of events before printout", (int*)&params->generation.gen_print_every );

      //-------------------------------------------------------------------------------------------
      // Process-specific parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<int>( "METH", "Computation method (kT-factorisation)", &proc_params_->operator[]<int>( "method" ) );
      registerParameter<int>( "IPOL", "Polarisation states to consider", &proc_params_->operator[]<int>( "polarisationStates" ) );

      //-------------------------------------------------------------------------------------------
      // Process kinematics parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<int>( "PMOD", "Outgoing primary particles' mode", &str_fun_ );
      registerParameter<int>( "EMOD", "Outgoing primary particles' mode", &str_fun_ );
      registerParameter<int>( "PAIR", "Outgoing particles' PDG id", (int*)&proc_params_->operator[]<int>( "pair" ) );

      registerParameter<int>( "INA1", "Heavy ion atomic weight (1st incoming beam)", (int*)&hi_1_.first );
      registerParameter<int>( "INZ1", "Heavy ion atomic number (1st incoming beam)", (int*)&hi_1_.second );
      registerParameter<int>( "INA2", "Heavy ion atomic weight (1st incoming beam)", (int*)&hi_2_.first );
      registerParameter<int>( "INZ2", "Heavy ion atomic number (1st incoming beam)", (int*)&hi_2_.second );
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
      try { setValue<double>( key.c_str(), std::stod( value ) ); } catch ( const std::invalid_argument& ) {}
      try { setValue<int>( key.c_str(), std::stoi( value ) ); } catch ( const std::invalid_argument& ) {}
      try { setValue<std::string>( key.c_str(), value ); } catch ( const std::invalid_argument& ) {
        throw CG_FATAL( "LpairHandler:setParameter" )
          << "Failed to add the parameter \"" << key << "\" → \"" << value << "\"!";
      }
    }

    std::string
    LpairHandler::getParameter( std::string key ) const
    {
      double dd = getValue<double>( key.c_str() );
      if ( dd != -999. )
        return std::to_string( dd );

      int ui = getValue<int>( key.c_str() );
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
