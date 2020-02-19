#include "CepGen/Cards/PythonHandler.h"

#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/EventModifier.h"

#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/ExportModule.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Processes/Process.h"
#include "CepGen/Modules/ProcessesFactory.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Integrator.h"

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/TamingFunction.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/PDG.h"

#include <algorithm>

#if PY_MAJOR_VERSION < 3
# define PYTHON2
#endif

namespace cepgen
{
  namespace card
  {
    PythonHandler::PythonHandler( const ParametersList& params ) :
      filename_( params.get<std::string>( FILENAME_KEY ) )
    {
      setenv( "PYTHONPATH", ".:Cards:test:../Cards", 1 );
      setenv( "PYTHONDONTWRITEBYTECODE", "1", 1 );
      CG_DEBUG( "PythonHandler" )
        << "Python PATH: " << getenv( "PYTHONPATH" ) << ".";
      if ( !filename_.empty() )
        parse( filename_ );
    }

    PythonHandler::PythonHandler( const std::string& file ) :
      filename_( file )
    {
      setenv( "PYTHONPATH", ".:Cards:test:../Cards", 1 );
      setenv( "PYTHONDONTWRITEBYTECODE", "1", 1 );
      CG_DEBUG( "PythonHandler" )
        << "Python PATH: " << getenv( "PYTHONPATH" ) << ".";
      if ( !filename_.empty() )
        parse( file );
    }

    Parameters&
    PythonHandler::parse( const std::string& file )
    {
      std::string filename = pythonPath( file );
      const size_t fn_len = filename.length()+1;

      //Py_DebugFlag = 1;
      //Py_VerboseFlag = 1;

      { // scope of the filename definition
#ifdef PYTHON2
        char* sfilename = new char[fn_len];
        snprintf( sfilename, fn_len, "%s", filename.c_str() );
#else
        wchar_t* sfilename = new wchar_t[fn_len];
        swprintf( sfilename, fn_len, L"%s", filename.c_str() );
#endif
        if ( !sfilename )
          throw CG_FATAL( "PythonHandler" )
            << "Invalid filename provided to the Python cards parser!";
        Py_SetProgramName( sfilename );
        delete [] sfilename;
      }

      Py_InitializeEx( 1 );

      if ( !Py_IsInitialized() )
        throw CG_FATAL( "PythonHandler" )
          << "Failed to initialise the Python cards parser!";

      CG_DEBUG( "PythonHandler" )
        << "Initialised the Python cards parser\n\t"
        << "Python version: " << Py_GetVersion() << "\n\t"
        << "Platform: " << Py_GetPlatform() << ".";

      PyObject* cfg = PyImport_ImportModule( filename.c_str() ); // new
      if ( !cfg )
        throwPythonError( "Failed to import the configuration card '"+file+"'" );

      //--- general particles definition
      if ( PyObject_HasAttrString( cfg, MCD_NAME ) == 1 ) {
        PyObject* ppdg = PyObject_GetAttrString( cfg, MCD_NAME ); // new
        if ( ppdg ) {
          pdg::MCDFileParser::parse( get<std::string>( ppdg ).c_str() );
          Py_CLEAR( ppdg );
        }
      }

      //--- additional particles definition
      if ( PyObject_HasAttrString( cfg, PDGLIST_NAME ) == 1 ) {
        PyObject* pextp = PyObject_GetAttrString( cfg, PDGLIST_NAME ); // new
        if ( pextp ) {
          parseExtraParticles( pextp );
          Py_CLEAR( pextp );
        }
      }

      //--- process definition
      PyObject* process = nullptr;
      if ( PyObject_HasAttrString( cfg, PROCESS_NAME ) != 1
        || !( process = PyObject_GetAttrString( cfg, PROCESS_NAME ) ) ) // new
        throwPythonError( "Failed to extract a '"+std::string( PROCESS_NAME )+"' keyword from the configuration card '"+file+"'!" );

      //--- list of process-specific parameters
      ParametersList proc_params;
      fillParameter( process, "processParameters", proc_params );

      //--- type of process to consider
      PyObject* pproc_name = element( process, ParametersList::MODULE_NAME ); // borrowed
      if ( !pproc_name )
        throwPythonError( "Failed to extract the process name from the configuration card '"+file+"'!" );
      const std::string proc_name = get<std::string>( pproc_name );

      //--- process mode
      params_.kinematics.mode = (KinematicsMode)proc_params.get<int>( "mode", (int)KinematicsMode::invalid );
      params_.setProcess( proc::ProcessesFactory::get().build( proc_name, proc_params ) );

      //--- process kinematics
      PyObject* pin_kinematics = element( process, "inKinematics" ); // borrowed
      if ( pin_kinematics )
        parseIncomingKinematics( pin_kinematics );

      PyObject* pout_kinematics = element( process, "outKinematics" ); // borrowed
      if ( pout_kinematics )
        parseOutgoingKinematics( pout_kinematics );

      //--- taming functions
      PyObject* ptam = element( process, "tamingFunctions" ); // borrowed
      if ( ptam )
        for ( const auto& p : getVector<ParametersList>( ptam ) )
          params_.taming_functions.emplace_back( p.get<std::string>( "variable" ), p.get<std::string>( "expression" ) );

      Py_CLEAR( process );

      if ( PyObject_HasAttrString( cfg, LOGGER_NAME ) == 1 ) {
        PyObject* plog = PyObject_GetAttrString( cfg, LOGGER_NAME ); // new
        if ( plog ) {
          parseLogging( plog );
          Py_CLEAR( plog );
        }
      }

      //--- hadroniser parameters (legacy)
      if ( PyObject_HasAttrString( cfg, HADR_NAME ) == 1 ) {
        PyObject* phad = PyObject_GetAttrString( cfg, HADR_NAME ); // new
        if ( phad ) {
          parseHadroniser( phad );
          Py_CLEAR( phad );
        }
      }

      if ( PyObject_HasAttrString( cfg, EVT_MOD_SEQ_NAME ) == 1 ) {
        PyObject* pmod_seq = PyObject_GetAttrString( cfg, EVT_MOD_SEQ_NAME ); // new
        if ( pmod_seq ) {
          parseEventModifiers( pmod_seq );
          Py_CLEAR( pmod_seq );
        }
      }

      //--- generation parameters
      if ( PyObject_HasAttrString( cfg, INTEGRATOR_NAME ) == 1 ) {
        PyObject* pint = PyObject_GetAttrString( cfg, INTEGRATOR_NAME ); // new
        if ( pint ) {
          parseIntegrator( pint );
          Py_CLEAR( pint );
        }
      }

      if ( PyObject_HasAttrString( cfg, GENERATOR_NAME ) == 1 ) {
        PyObject* pgen = PyObject_GetAttrString( cfg, GENERATOR_NAME ); // new
        if ( pgen ) {
          parseGenerator( pgen );
          Py_CLEAR( pgen );
        }
      }

      if ( PyObject_HasAttrString( cfg, OUTPUT_NAME ) == 1 ) {
        PyObject* pout = PyObject_GetAttrString( cfg, OUTPUT_NAME ); // new
        if ( pout ) {
          if ( isVector<ParametersList>( pout ) )
            parseOutputModules( pout );
          else
            parseOutputModule( pout );
          Py_CLEAR( pout );
        }
      }

      //--- finalisation
      Py_CLEAR( cfg );

      if ( Py_IsInitialized() )
        Py_Finalize();

      return params_;
    }

    void
    PythonHandler::parseIncomingKinematics( PyObject* kin )
    {
      //--- retrieve the beams PDG ids
      std::vector<ParametersList> beams_pdg;
      fillParameter( kin, "pdgIds", beams_pdg );
      if ( !beams_pdg.empty() ) {
        if ( beams_pdg.size() != 2 )
          throwPythonError( utils::format( "Invalid list of PDG ids retrieved for incoming beams:\n\t2 PDG ids are expected, %d provided!", beams_pdg.size() ) );
        params_.kinematics.incoming_beams. first.pdg = (pdgid_t)beams_pdg.at( 0 ).get<int>( "pdgid" );
        params_.kinematics.incoming_beams.second.pdg = (pdgid_t)beams_pdg.at( 1 ).get<int>( "pdgid" );
      }
      //--- incoming beams kinematics
      std::vector<double> beams_pz;
      fillParameter( kin, "pz", beams_pz );
      if ( !beams_pz.empty() ) {
        if ( beams_pz.size() != 2 )
          throwPythonError( utils::format( "Invalid list of pz's retrieved for incoming beams:\n\t2 pz's are expected, %d provided!", beams_pz.size() ) );
        params_.kinematics.incoming_beams. first.pz = beams_pz.at( 0 );
        params_.kinematics.incoming_beams.second.pz = beams_pz.at( 1 );
      }
      double sqrt_s = -1.;
      fillParameter( kin, "cmEnergy", sqrt_s );
      if ( sqrt_s != -1. )
        params_.kinematics.setSqrtS( sqrt_s );
      //--- structure functions set for incoming beams
      PyObject* psf = element( kin, "structureFunctions" ); // borrowed
      if ( psf )
        params_.kinematics.structure_functions = strfun::StructureFunctionsFactory::get().build( get<ParametersList>( psf ) );
      //--- types of parton fluxes for kt-factorisation
      PyObject* pktf = element( kin, "ktFluxes" ); // borrowed
      if ( pktf ) {
        if ( isVector<int>( pktf ) ) {
          std::vector<int> kt_fluxes;
          fillParameter( kin, "ktFluxes", kt_fluxes );
          if ( !kt_fluxes.empty() ) {
            params_.kinematics.incoming_beams.first.kt_flux = (KTFlux)kt_fluxes.at( 0 );
            params_.kinematics.incoming_beams.second.kt_flux = ( kt_fluxes.size() > 1 )
              ? (KTFlux)kt_fluxes.at( 1 )
              : (KTFlux)kt_fluxes.at( 0 );
          }
        }
        else if ( is<int>( pktf ) ) {
          int kt_fluxes;
          fillParameter( kin, "ktFluxes", kt_fluxes );
          params_.kinematics.incoming_beams.first.kt_flux = params_.kinematics.incoming_beams.second.kt_flux = (KTFlux)kt_fluxes;
        }
        else
          throwPythonError( "Unsupported format for the ktFluxes definition!" );
      }
      //--- specify where to look for the grid path for gluon emission
      std::string kmr_grid_path;
      fillParameter( kin, "kmrGridPath", kmr_grid_path );
      if ( !kmr_grid_path.empty() )
        kmr::GluonGrid::get( kmr_grid_path.c_str() );
      //--- parse heavy ions beams
      std::vector<int> hi_beam1, hi_beam2;
      fillParameter( kin, "heavyIonA", hi_beam1 );
      if ( hi_beam1.size() == 2 )
        params_.kinematics.incoming_beams. first.pdg = HeavyIon{ (unsigned short)hi_beam1[0], (Element)hi_beam1[1] };
      fillParameter( kin, "heavyIonB", hi_beam2 );
      if ( hi_beam2.size() == 2 )
        params_.kinematics.incoming_beams.second.pdg = HeavyIon{ (unsigned short)hi_beam2[0], (Element)hi_beam2[1] };
    }

    void
    PythonHandler::parseOutgoingKinematics( PyObject* kin )
    {
      std::vector<int> parts;
      fillParameter( kin, "minFinalState", parts );
      for ( const auto& pdg : parts )
        params_.kinematics.minimum_final_state.emplace_back( (pdgid_t)pdg );

      ParametersList part_cuts;
      fillParameter( kin, "cuts", part_cuts );
      for ( const auto& part : part_cuts.keys() ) {
        const auto pdg = (pdgid_t)stoi( part );
        const auto& cuts = part_cuts.get<ParametersList>( part );
        if ( cuts.has<Limits>( "pt" ) )
          params_.kinematics.cuts.central_particles[pdg].pt_single = cuts.get<Limits>( "pt" );
        if ( cuts.has<Limits>( "energy" ) )
          params_.kinematics.cuts.central_particles[pdg].energy_single = cuts.get<Limits>( "energy" );
        if ( cuts.has<Limits>( "eta" ) )
          params_.kinematics.cuts.central_particles[pdg].eta_single = cuts.get<Limits>( "eta" );
        if ( cuts.has<Limits>( "rapidity" ) )
          params_.kinematics.cuts.central_particles[pdg].rapidity_single = cuts.get<Limits>( "rapidity" );
      }

      // for LPAIR/collinear matrix elements
      fillParameter( kin, "q2", params_.kinematics.cuts.initial.q2 );

      // for the kT factorised matrix elements
      fillParameter( kin, "qt", params_.kinematics.cuts.initial.qt );
      fillParameter( kin, "phiqt", params_.kinematics.cuts.initial.phi_qt );
      fillParameter( kin, "ptdiff", params_.kinematics.cuts.central.pt_diff );
      fillParameter( kin, "phiptdiff", params_.kinematics.cuts.central.phi_pt_diff );
      fillParameter( kin, "rapiditydiff", params_.kinematics.cuts.central.rapidity_diff );

      // generic phase space limits
      fillParameter( kin, "rapidity", params_.kinematics.cuts.central.rapidity_single );
      fillParameter( kin, "eta", params_.kinematics.cuts.central.eta_single );
      fillParameter( kin, "pt", params_.kinematics.cuts.central.pt_single );

      fillParameter( kin, "ptsum", params_.kinematics.cuts.central.pt_sum );
      fillParameter( kin, "invmass", params_.kinematics.cuts.central.mass_sum );

      fillParameter( kin, "mx", params_.kinematics.cuts.remnants.mass_single );
      fillParameter( kin, "yj", params_.kinematics.cuts.remnants.rapidity_single );

      Limits lim_xi;
      fillParameter( kin, "xi", lim_xi );
      if ( lim_xi.valid() )
        //params_.kinematics.cuts.remnants.energy_single = ( lim_xi+(-1.) )*( -params_.kinematics.incoming_beams.first.pz );
        params_.kinematics.cuts.remnants.energy_single = -( lim_xi-1. )*params_.kinematics.incoming_beams.first.pz;
    }

    void
    PythonHandler::parseLogging( PyObject* log )
    {
      int log_level = 0;
      fillParameter( log, "level", log_level );
      utils::Logger::get().level = (utils::Logger::Level)log_level;
      std::vector<std::string> enabled_modules;
      fillParameter( log, "enabledModules", enabled_modules );
      for ( const auto& mod : enabled_modules )
        utils::Logger::get().addExceptionRule( mod );
    }

    void
    PythonHandler::parseIntegrator( PyObject* integr )
    {
      if ( !PyDict_Check( integr ) )
        throwPythonError( "Integrator object should be a dictionary!" );
      PyObject* palgo = element( integr, ParametersList::MODULE_NAME ); // borrowed
      if ( !palgo )
        throwPythonError( "Failed to retrieve the integration algorithm name!" );
      std::string algo = get<std::string>( palgo );
      if ( algo == "plain" )
        params_.integration().type = IntegratorType::plain;
      else if ( algo == "Vegas" ) {
        params_.integration().type = IntegratorType::Vegas;
        fillParameter( integr, "alpha", (double&)params_.integration().vegas.alpha );
        fillParameter( integr, "iterations", params_.integration().vegas.iterations );
        fillParameter( integr, "mode", (int&)params_.integration().vegas.mode );
        fillParameter( integr, "verbosity", (int&)params_.integration().vegas.verbose );
        std::string vegas_logging_output = "cerr";
        fillParameter( integr, "loggingOutput", vegas_logging_output );
        if ( vegas_logging_output == "cerr" )
          // redirect all debugging information to the error stream
          params_.integration().vegas.ostream = stderr;
        else if ( vegas_logging_output == "cout" )
          // redirect all debugging information to the standard stream
          params_.integration().vegas.ostream = stdout;
        else
          params_.integration().vegas.ostream = fopen( vegas_logging_output.c_str(), "w" );
      }
      else if ( algo == "MISER" ) {
        params_.integration().type = IntegratorType::MISER;
        fillParameter( integr, "estimateFraction", (double&)params_.integration().miser.estimate_frac );
        fillParameter( integr, "minCalls", params_.integration().miser.min_calls );
        fillParameter( integr, "minCallsPerBisection", params_.integration().miser.min_calls_per_bisection );
        fillParameter( integr, "alpha", (double&)params_.integration().miser.alpha );
        fillParameter( integr, "dither", (double&)params_.integration().miser.dither );
      }
      else
        throwPythonError( utils::format( "Invalid integration() algorithm: %s", algo.c_str() ) );

      fillParameter( integr, "numFunctionCalls", params_.integration().ncvg );
      fillParameter( integr, "seed", (unsigned long&)params_.integration().rng_seed );
      unsigned int rng_engine;
      fillParameter( integr, "rngEngine", rng_engine );
      switch ( rng_engine ) {
        case 0: default: params_.integration().rng_engine = (gsl_rng_type*)gsl_rng_mt19937; break;
        case 1: params_.integration().rng_engine = (gsl_rng_type*)gsl_rng_taus2; break;
        case 2: params_.integration().rng_engine = (gsl_rng_type*)gsl_rng_gfsr4; break;
        case 3: params_.integration().rng_engine = (gsl_rng_type*)gsl_rng_ranlxs0; break;
      }
      fillParameter( integr, "chiSqCut", params_.integration().vegas_chisq_cut );
    }

    void
    PythonHandler::parseGenerator( PyObject* gen )
    {
      if ( !PyDict_Check( gen ) )
        throwPythonError( "Generation information object should be a dictionary!" );
      params_.generation().enabled = true;
      fillParameter( gen, "treat", params_.generation().treat );
      fillParameter( gen, "numEvents", params_.generation().maxgen );
      fillParameter( gen, "printEvery", params_.generation().gen_print_every );
      fillParameter( gen, "numThreads", params_.generation().num_threads );
      fillParameter( gen, "numPoints", params_.generation().num_points );
    }

    void
    PythonHandler::parseEventModifiers( PyObject* mod )
    {
      if ( !PyList_Check( mod ) )
        throwPythonError( "Event modification definition object should be a list/Sequence!" );

      for ( Py_ssize_t i = 0; i < PyList_Size( mod ); ++i )
        parseHadroniser( PyList_GetItem( mod, i ) );
    }

    void
    PythonHandler::parseHadroniser( PyObject* mod )
    {
      if ( !PyDict_Check( mod ) )
        throwPythonError( "Event modification definition object should be a dictionary!" );

      PyObject* pname = element( mod, ParametersList::MODULE_NAME ); // borrowed
      if ( !pname )
        throwPythonError( "Event modification algorithm name is required!" );
      std::string mod_name = get<std::string>( pname );

      params_.addModifier( EventModifierFactory::get().build( mod_name, get<ParametersList>( mod ) ) );

      auto h = params_.eventModifiersSequence().rbegin()->get();
      h->setParameters( params_ );
      { //--- before calling the init() method
        std::vector<std::string> config;
        fillParameter( mod, "preConfiguration", config );
        h->readStrings( config );
      }
      h->init();
      { //--- after init() has been called
        std::vector<std::string> config;
        fillParameter( mod, "processConfiguration", config );
        for ( const auto& block : config ) {
          std::vector<std::string> config_blk;
          fillParameter( mod, block.c_str(), config_blk );
          h->readStrings( config_blk );
        }
      }
    }

    void
    PythonHandler::parseOutputModules( PyObject* mod )
    {
      if ( !PyList_Check( mod ) )
        throwPythonError( "Output modules definition object should be a list/Sequence!" );

      for ( Py_ssize_t i = 0; i < PyList_Size( mod ); ++i )
        parseOutputModule( PyList_GetItem( mod, i ) );
    }

    void
    PythonHandler::parseOutputModule( PyObject* pout )
    {
      if ( !is<ParametersList>( pout ) )
        throwPythonError( "Invalid type for output parameters list!" );

      PyObject* pname = element( pout, ParametersList::MODULE_NAME ); // borrowed
      if ( !pname )
        throwPythonError( "Output module name is required!" );
      params_.addOutputModule( io::ExportModuleFactory::get().build( get<std::string>( pname ), get<ParametersList>( pout ) ) );
    }

    void
    PythonHandler::parseExtraParticles( PyObject* pparts )
    {
      if ( !is<ParametersList>( pparts ) )
        throwPythonError( "Extra particles definition object should be a parameters list!" );

      const auto& parts = get<ParametersList>( pparts );
      for ( const auto& k : parts.keys() ) {
        const auto& part = parts.get<ParticleProperties>( k );
        if ( part.pdgid == 0 || part.mass < 0. )
          continue;
        CG_DEBUG( "PythonHandler:particles" )
          << "Adding a new particle with name \"" << part.name << "\" to the PDG dictionary.";
        PDG::get().define( part );
      }
    }
  }
}

REGISTER_CARD_HANDLER( "py", PythonHandler )

