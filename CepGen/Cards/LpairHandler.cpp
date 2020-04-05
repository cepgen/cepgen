#include "CepGen/Cards/LpairHandler.h"

#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Modules/EventModifierFactory.h"

#include "CepGen/Core/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Modules/ProcessesFactory.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/PDG.h"

#include <fstream>
#include <iomanip>

namespace cepgen
{
  namespace card
  {
    const int LpairHandler::kInvalid = 99999;

    //----- specialization for LPAIR input cards

    LpairHandler::LpairHandler( const ParametersList& params ) :
      proc_params_( new ParametersList ),
      str_fun_( 11 ), sr_type_( 1 ), xi_min_( 0. ), xi_max_( 1. ),
      pdg_input_path_( "External/mass_width_2019.mcd" ), iend_( 1 ),
      hi_1_( { 0, 0 } ), hi_2_( { 0, 0 } )
    {
      const auto file = params.get<std::string>( FILENAME_KEY );
      if ( !file.empty() )
        parse( file, params_ );
    }

    Parameters
    LpairHandler::parse( const std::string& file, Parameters& params )
    {
      std::ifstream f( file, std::fstream::in );
      if ( !f.is_open() )
        throw CG_FATAL( "LpairHandler" ) << "Failed to parse file \"" << file << "%s\".";

      init();

      //--- parse all fields
      std::string key, value;
      std::ostringstream os;
      while ( f >> key >> value ) {
        if ( key[0] == '#' ) // FIXME need to ensure there is no extra space before!
          continue;
        setParameter( key, value );
        if ( description( key ) != "null" )
          os << "\n>> " << key << " = " << std::setw( 25 ) << parameter( key )
             << " (" << description( key ) << ")";
      }
      f.close();
      params_ = params;

      CG_INFO( "LpairHandler" ) << "File '" << file << "' succesfully opened!\n\t"
        << "The following parameters are set:" << os.str() << "\n\t"
        << "Now parsing the configuration.";

      //--- parse the PDG library
      if ( !pdg_input_path_.empty() )
        pdg::MCDFileParser::parse( pdg_input_path_.c_str() );
      if ( !kmr_grid_path_.empty() )
        kmr::GluonGrid::get( kmr_grid_path_.c_str() );

      //--- parse the process name
      if ( !proc_name_.empty() || !proc_params_->empty() ) {
        if ( !params.hasProcess() && proc_name_.empty() )
          throw CG_FATAL( "LpairHandler" ) << "Process name not specified!";
        if ( params_.hasProcess() && params_.process().name() == proc_name_ )
          *proc_params_ = ParametersList( params_.process().parameters() )+*proc_params_;
        params_.setProcess( proc::ProcessesFactory::get().build( proc_name_, *proc_params_ ) );
      }

      const Limits lim_xi{ xi_min_, xi_max_ };
      if ( lim_xi.valid() )
        params_.kinematics.cuts.remnants.energy_single = ( lim_xi+(-1.) )*( -params_.kinematics.incoming_beams.first.pz );

      //--- parse the structure functions code
      auto sf_params = ParametersList()
        .setName<int>( str_fun_ )
        .set<ParametersList>( "sigmaRatio", ParametersList()
          .setName<int>( sr_type_ ) );
      const unsigned long kLHAPDFCodeDec = 10000000, kLHAPDFPartDec = 1000000;
      if ( str_fun_ / kLHAPDFCodeDec == 1 ) { // SF from parton
        const unsigned long icode = str_fun_ % kLHAPDFCodeDec;
        sf_params
          .setName<int>( (int)strfun::Type::Partonic )
          .set<int>( "pdfId", icode % kLHAPDFPartDec )
          .set<int>( "mode", icode / kLHAPDFPartDec ); // 0, 1, 2
      }
      else if ( str_fun_ == (int)strfun::Type::MSTWgrid )
        sf_params
          .set<std::string>( "gridPath", mstw_grid_path_ );
      params_.kinematics.structure_functions = strfun::StructureFunctionsFactory::get().build( sf_params );

      //--- check if event generation is required
      params_.generation().enabled = iend_ > 1;

      //--- parse the hadronisation algorithm name
      if ( !evt_mod_name_.empty() )
        for ( const auto& mod : utils::split( evt_mod_name_, ',' ) ) {
          params_.addModifier( EventModifierFactory::get().build( mod, ParametersList() ) );
          (*params_.eventModifiersSequence().rbegin())->setParameters( params_ );
        }

      //--- parse the output module name
      if ( !out_mod_name_.empty() ) {
        const auto& out_files = utils::split( out_file_name_, ',' );
        size_t i = 0;
        for ( const auto& mod : utils::split( out_mod_name_, ',' ) ) {
          ParametersList outm;
          if ( out_files.size() > i && !out_files.at( i ).empty() )
            outm.set<std::string>( "filename", out_files.at( i ) );
          params_.addOutputModule( io::ExportModuleFactory::get().build( mod, outm ) );
          ++i;
        }
      }

      //--- check if we are dealing with heavy ions for incoming states
      HeavyIon hi1{ hi_1_.first, (Element)hi_1_.second }, hi2{ hi_2_.first, (Element)hi_2_.second };
      if ( hi1 )
        params_.kinematics.incoming_beams.first.pdg = hi1;
      if ( hi2 )
        params_.kinematics.incoming_beams.second.pdg = hi2;

      return params_;
    }

    void
    LpairHandler::init()
    {
      //-------------------------------------------------------------------------------------------
      // Process/integration/hadronisation parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<std::string>( "PROC", "Process name to simulate", &proc_name_ );
      registerParameter<std::string>( "ITYP", "Integration algorithm", &params_.integrator->operator[]<std::string>( ParametersList::MODULE_NAME ) );
      registerParameter<std::string>( "HADR", "Hadronisation algorithm", &evt_mod_name_ );
      registerParameter<std::string>( "EVMD", "Events modification algorithms", &evt_mod_name_ );
      registerParameter<std::string>( "OUTP", "Output module", &out_mod_name_ );
      registerParameter<std::string>( "OUTF", "Output file name", &out_file_name_ );

      //-------------------------------------------------------------------------------------------
      // General parameters
      //-------------------------------------------------------------------------------------------

      registerParameter<bool>( "NTRT", "Smoothen the integrand", &params_.integrator->operator[]<bool>( "treat" ) );
      registerParameter<int>( "IEND", "Generation type", &iend_ );
      registerParameter<int>( "DEBG", "Debugging verbosity", (int*)&utils::Logger::get().level );
      registerParameter<int>( "NCVG", "Number of function calls", (int*)&params_.integrator->operator[]<int>( "numFunctionCalls" ) );
      registerParameter<int>( "ITVG", "Number of integration iterations", (int*)&params_.integrator->operator[]<int>( "iterations" ) );
      registerParameter<int>( "SEED", "Random generator seed", (int*)&params_.integrator->operator[]<int>( "seed" ) );
      registerParameter<int>( "NTHR", "Number of threads to use for events generation", (int*)&params_.generation().num_threads );
      registerParameter<int>( "MODE", "Subprocess' mode", (int*)&params_.kinematics.mode );
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
      registerParameter<std::string>( "PDGI", "Input file for PDG information", &pdg_input_path_ );
      registerParameter<int>( "PMOD", "Outgoing primary particles' mode", &str_fun_ );
      registerParameter<int>( "EMOD", "Outgoing primary particles' mode", &str_fun_ );
      registerParameter<int>( "RTYP", "R-ratio computation type", &sr_type_ );
      registerParameter<int>( "PAIR", "Outgoing particles' PDG id", (int*)&proc_params_->operator[]<int>( "pair" ) );
      registerParameter<int>( "INA1", "Heavy ion atomic weight (1st incoming beam)", (int*)&hi_1_.first );
      registerParameter<int>( "INZ1", "Heavy ion atomic number (1st incoming beam)", (int*)&hi_1_.second );
      registerParameter<int>( "INA2", "Heavy ion atomic weight (2nd incoming beam)", (int*)&hi_2_.first );
      registerParameter<int>( "INZ2", "Heavy ion atomic number (2nd incoming beam)", (int*)&hi_2_.second );
      registerParameter<double>( "INP1", "Momentum (1st primary particle)", &params_.kinematics.incoming_beams.first.pz );
      registerParameter<double>( "INP2", "Momentum (2nd primary particle)", &params_.kinematics.incoming_beams.second.pz );
      registerParameter<double>( "INPP", "Momentum (1st primary particle)", &params_.kinematics.incoming_beams.first.pz );
      registerParameter<double>( "INPE", "Momentum (2nd primary particle)", &params_.kinematics.incoming_beams.second.pz );
      registerParameter<double>( "PTCT", "Minimal transverse momentum (single central outgoing particle)", &params_.kinematics.cuts.central.pt_single.min() );
      registerParameter<double>( "PTMX", "Maximal transverse momentum (single central outgoing particle)", &params_.kinematics.cuts.central.pt_single.max() );
      registerParameter<double>( "MSCT", "Minimal central system mass", &params_.kinematics.cuts.central.mass_sum.min() );
      registerParameter<double>( "ECUT", "Minimal energy (single central outgoing particle)", &params_.kinematics.cuts.central.energy_single.min() );
      registerParameter<double>( "ETMN", "Minimal pseudo-rapidity (central outgoing particles)", &params_.kinematics.cuts.central.eta_single.min() );
      registerParameter<double>( "ETMX", "Maximal pseudo-rapidity (central outgoing particles)", &params_.kinematics.cuts.central.eta_single.max() );
      registerParameter<double>( "YMIN", "Minimal rapidity (central outgoing particles)", &params_.kinematics.cuts.central.rapidity_single.min() );
      registerParameter<double>( "YMAX", "Maximal rapidity (central outgoing particles)", &params_.kinematics.cuts.central.rapidity_single.max() );
      registerParameter<double>( "PDMN", "Minimal transverse momentum difference (central outgoing particles)", &params_.kinematics.cuts.central.pt_diff.min() );
      registerParameter<double>( "PDMX", "Maximal transverse momentum difference (central outgoing particles)", &params_.kinematics.cuts.central.pt_diff.max() );
      registerParameter<double>( "Q2MN", "Minimal Q² = -q² (exchanged parton)", &params_.kinematics.cuts.initial.q2.min() );
      registerParameter<double>( "Q2MX", "Maximal Q² = -q² (exchanged parton)", &params_.kinematics.cuts.initial.q2.max() );
      registerParameter<double>( "QTMN", "Minimal Q_T (exchanged parton)", &params_.kinematics.cuts.initial.qt.min() );
      registerParameter<double>( "QTMX", "Maximal Q_T (exchanged parton)", &params_.kinematics.cuts.initial.qt.max() );
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
      f.close();
    }

    void
    LpairHandler::setParameter( const std::string& key, const std::string& value )
    {
      // particular case for the double as we cannot rely on casting exceptions
      if ( value.find( '.' ) != std::string::npos )
        try {
          setValue<double>( key.c_str(), std::stod( value ) );
          return;
        } catch ( const std::logic_error& ) {
          throw CG_FATAL( "LpairHandler:setParameter" )
            << "Failed to parse a floating-point parameter \"" << key << "\" → \"" << value << "\"!";
        }
      try { setValue<int>( key.c_str(), std::stoi( value ) ); } catch ( const std::logic_error& ) {
        try { setValue<std::string>( key.c_str(), value ); } catch ( const std::logic_error& ) {
          throw CG_FATAL( "LpairHandler:setParameter" )
            << "Failed to add the parameter \"" << key << "\" → \"" << value << "\"!";
        }
      }
    }

    std::string
    LpairHandler::parameter( std::string key ) const
    {
      {
        auto var = getValue<double>( key.c_str() );
        if ( var != -999. )
          return std::to_string( var );
      }{
        auto var = getValue<int>( key.c_str() );
        if ( var != -999999 )
          return std::to_string( var );
      }
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
      return "null";
    }
  }
}

REGISTER_CARD_HANDLER( "card", LpairHandler )
