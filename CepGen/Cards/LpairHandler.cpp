#include "CepGen/Cards/LpairHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Integrator.h"

#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/Core/EventModifierHandler.h"
#include "CepGen/IO/ExportHandler.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/FormFactors.h"

#include <fstream>
#include <iomanip>

namespace cepgen
{
  namespace card
  {
    const int LpairHandler::kInvalid = 99999;

    //----- specialization for LPAIR input cards

    LpairHandler::LpairHandler( const char* file ) :
      proc_params_( new ParametersList ),
      str_fun_( 11 ), sr_type_( 1 ), xi_min_( 0. ), xi_max_( 1. ),
      hi_1_( { 0, 0 } ), hi_2_( { 0, 0 } )
    {
      std::ifstream f( file, std::fstream::in );
      if ( !f.is_open() )
        throw CG_FATAL( "LpairHandler" ) << "Failed to parse file \"" << file << "%s\".";

      init();

      //--- parse all fields
      std::unordered_map<std::string, std::string> m_params;
      std::string key, value;
      std::ostringstream os;
      while ( f >> key >> value ) {
        if ( key[0] == '#' ) // FIXME need to ensure there is no extra space before!
          continue;
        setParameter( key, value );
        m_params.insert( { key, value } );
        if ( description( key ) != "null" )
          os << "\n>> " << key << " = " << std::setw( 15 ) << parameter( key )
             << " (" << description( key ) << ")";
      }
      f.close();

      //--- parse the process name
      params_.setProcess( proc::ProcessesHandler::get().build( proc_name_, *proc_params_ ) );

      const Limits lim_xi{ xi_min_, xi_max_ };
      if ( lim_xi.valid() )
        params_.kinematics.cuts.remnants.energy_single = ( lim_xi+(-1.) )*( -params_.kinematics.incoming_beams.first.pz );

      //--- parse the structure functions code
      auto sf_params = ParametersList()
        .set<int>( strfun::StructureFunctionsHandler::KEY, str_fun_ )
        .set<ParametersList>( "sigmaRatio", ParametersList()
          .set<int>( "id", sr_type_ ) );
      const unsigned long kLHAPDFCodeDec = 10000000, kLHAPDFPartDec = 1000000;
      if ( str_fun_ / kLHAPDFCodeDec == 1 ) { // SF from parton
        const unsigned long icode = str_fun_ % kLHAPDFCodeDec;
        sf_params
          .set<int>( strfun::StructureFunctionsHandler::KEY, (int)strfun::Type::Partonic )
          .set<int>( "pdfId", icode % kLHAPDFPartDec )
          .set<int>( "mode", icode / kLHAPDFPartDec ); // 0, 1, 2
      }
      else if ( str_fun_ == (int)strfun::Type::MSTWgrid )
        sf_params
          .set<std::string>( "gridPath", mstw_grid_path_ );
      params_.kinematics.incoming_beams.first.form_factors->setStructureFunctions( strfun::StructureFunctionsHandler::get().build( sf_params ) );
      params_.kinematics.incoming_beams.second.form_factors->setStructureFunctions( strfun::StructureFunctionsHandler::get().build( sf_params ) );

      //--- parse the integration algorithm name
      if ( integr_type_ == "plain" )
        params_.integration().type = IntegratorType::plain;
      else if ( integr_type_ == "Vegas" )
        params_.integration().type = IntegratorType::Vegas;
      else if ( integr_type_ == "MISER" )
        params_.integration().type = IntegratorType::MISER;
      else if ( integr_type_ != "" )
        throw CG_FATAL( "LpairHandler" ) << "Unrecognized integrator type: " << integr_type_ << "!";

      //--- parse the hadronisation algorithm name
      if ( !evt_mod_name_.empty() )
        for ( const auto& mod : split( evt_mod_name_, ',' ) ) {
          params_.addModifier( cepgen::EventModifierHandler::get().build( mod, ParametersList() ) );
          (*params_.eventModifiersSequence().rbegin())->setParameters( params_ );
        }

      //--- parse the output module name
      if ( !out_mod_name_.empty() ) {
        ParametersList outm;
        if ( !out_file_name_.empty() )
          outm.set<std::string>( "filename", out_file_name_ );
        params_.setOutputModule( cepgen::io::ExportHandler::get().build( out_mod_name_, outm ) );
      }

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
    LpairHandler::init()
    {
      //-------------------------------------------------------------------------------------------
      // Process/integration/hadronisation parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>( "PROC", "Process name to simulate", &proc_name_ );
      registerParameter<std::string>( "ITYP", "Integration algorithm", &integr_type_ );
      registerParameter<std::string>( "HADR", "Hadronisation algorithm", &evt_mod_name_ );
      registerParameter<std::string>( "EVMD", "Events modification algorithms", &evt_mod_name_ );
      registerParameter<std::string>( "OUTP", "Output module", &out_mod_name_ );
      registerParameter<std::string>( "OUTF", "Output file name", &out_file_name_ );

      //-------------------------------------------------------------------------------------------
      // General parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<bool>( "IEND", "Generation type", &params_.generation().enabled );
      registerParameter<bool>( "NTRT", "Smoothen the integrand", &params_.generation().treat );
      registerParameter<int>( "DEBG", "Debugging verbosity", (int*)&utils::Logger::get().level );
      registerParameter<int>( "NCVG", "Number of function calls", (int*)&params_.integration().ncvg );
      registerParameter<int>( "ITVG", "Number of integration iterations", (int*)&params_.integration().vegas.iterations );
      registerParameter<int>( "SEED", "Random generator seed", (int*)&params_.integration().rng_seed );
      registerParameter<int>( "NTHR", "Number of threads to use for events generation", (int*)&params_.generation().num_threads );
      registerParameter<int>( "MODE", "Subprocess' mode", (int*)&proc_params_->operator[]<int>( "mode" ) );
      registerParameter<int>( "NCSG", "Number of points to probe", (int*)&params_.generation().num_points );
      registerParameter<int>( "NGEN", "Number of events to generate", (int*)&params_.generation().maxgen );
      registerParameter<int>( "NPRN", "Number of events before printout", (int*)&params_.generation().gen_print_every );

      //-------------------------------------------------------------------------------------------
      // Process-specific parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<int>( "METH", "Computation method (kT-factorisation)", &proc_params_->operator[]<int>( "method" ) );
      registerParameter<int>( "IPOL", "Polarisation states to consider", &proc_params_->operator[]<int>( "polarisationStates" ) );

      //-------------------------------------------------------------------------------------------
      // Process kinematics parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>( "KMRG", "KMR grid interpolation path", &kmr_grid_path_ );
      registerParameter<std::string>( "MGRD", "MSTW grid interpolation path", &mstw_grid_path_ );
      registerParameter<int>( "PMOD", "Outgoing primary particles' mode", &str_fun_ );
      registerParameter<int>( "EMOD", "Outgoing primary particles' mode", &str_fun_ );
      registerParameter<int>( "RTYP", "R-ratio computation type", &sr_type_ );
      registerParameter<int>( "PAIR", "Outgoing particles' PDG id", (int*)&proc_params_->operator[]<int>( "pair" ) );
      registerParameter<int>( "INA1", "Heavy ion atomic weight (1st incoming beam)", (int*)&hi_1_.first );
      registerParameter<int>( "INZ1", "Heavy ion atomic number (1st incoming beam)", (int*)&hi_1_.second );
      registerParameter<int>( "INA2", "Heavy ion atomic weight (1st incoming beam)", (int*)&hi_2_.first );
      registerParameter<int>( "INZ2", "Heavy ion atomic number (1st incoming beam)", (int*)&hi_2_.second );
      registerParameter<double>( "INP1", "Momentum (1st primary particle)", &params_.kinematics.incoming_beams.first.pz );
      registerParameter<double>( "INP2", "Momentum (2nd primary particle)", &params_.kinematics.incoming_beams.second.pz );
      registerParameter<double>( "INPP", "Momentum (1st primary particle)", &params_.kinematics.incoming_beams.first.pz );
      registerParameter<double>( "INPE", "Momentum (2nd primary particle)", &params_.kinematics.incoming_beams.second.pz );
      registerParameter<double>( "PTCT", "Minimal transverse momentum (single central outgoing particle)", &params_.kinematics.cuts.central.pt_single.min() );
      registerParameter<double>( "MSCT", "Minimal central system mass", &params_.kinematics.cuts.central.mass_sum.min() );
      registerParameter<double>( "ECUT", "Minimal energy (single central outgoing particle)", &params_.kinematics.cuts.central.energy_single.min() );
      registerParameter<double>( "ETMN", "Minimal pseudo-rapidity (central outgoing particles)", &params_.kinematics.cuts.central.eta_single.min() );
      registerParameter<double>( "ETMX", "Maximal pseudo-rapidity (central outgoing particles)", &params_.kinematics.cuts.central.eta_single.max() );
      registerParameter<double>( "YMIN", "Minimal rapidity (central outgoing particles)", &params_.kinematics.cuts.central.rapidity_single.min() );
      registerParameter<double>( "YMAX", "Maximal rapidity (central outgoing particles)", &params_.kinematics.cuts.central.rapidity_single.max() );
      registerParameter<double>( "Q2MN", "Minimal Q² = -q² (exchanged parton)", &params_.kinematics.cuts.initial.q2.min() );
      registerParameter<double>( "Q2MX", "Maximal Q² = -q² (exchanged parton)", &params_.kinematics.cuts.initial.q2.max() );
      registerParameter<double>( "MXMN", "Minimal invariant mass of proton remnants", &params_.kinematics.cuts.remnants.mass_single.min() );
      registerParameter<double>( "MXMX", "Maximal invariant mass of proton remnants", &params_.kinematics.cuts.remnants.mass_single.max() );
      registerParameter<double>( "XIMN", "Minimal fractional momentum loss of outgoing proton (ξ)", &xi_min_ );
      registerParameter<double>( "XIMX", "Maximal fractional momentum loss of outgoing proton (ξ)", &xi_max_ );
    }

    void
    LpairHandler::store( const char* file )
    {
      std::ofstream f( file, std::fstream::out | std::fstream::trunc );
      if ( !f.is_open() )
        throw CG_ERROR( "LpairHandler" ) << "Failed to open file \"" << file << "%s\" for writing.";
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
      try { setValue<double>( key.c_str(), std::stod( value ) ); } catch ( const std::invalid_argument& ) {
        try { setValue<int>( key.c_str(), std::stoi( value ) ); } catch ( const std::invalid_argument& ) {
          try { setValue<std::string>( key.c_str(), value ); } catch ( const std::invalid_argument& ) {
            throw CG_FATAL( "LpairHandler:setParameter" )
              << "Failed to add the parameter \"" << key << "\" → \"" << value << "\"!";
          }
        }
      }
    }

    std::string
    LpairHandler::parameter( std::string key ) const
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
    LpairHandler::description( std::string key ) const
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

    std::vector<std::string>
    LpairHandler::split( const std::string& str, char delim )
    {
      std::vector<std::string> out;
      std::string token;
      std::istringstream iss( str );
      while ( std::getline( iss, token, delim ) )
        out.emplace_back( token );
      return out;
    }
  }
}
