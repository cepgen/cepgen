#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"

#ifdef PYTHON

#include "CepGen/Core/TamingFunction.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoFF.h"
#include "CepGen/Processes/PPtoWW.h"
#include "CepGen/Processes/FortranProcesses.h"

#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/StructureFunctions/LHAPDF.h"
#include "CepGen/StructureFunctions/MSTWGrid.h"
#include "CepGen/StructureFunctions/Schaefer.h"

#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

#include <algorithm>

#if PY_MAJOR_VERSION < 3
#  define PYTHON2
#endif

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for CepGen input cards
    PythonHandler::PythonHandler( const char* file )
    {
      setenv( "PYTHONPATH", ".:..:Cards", 1 );
      std::string filename = getPythonPath( file );
      const size_t fn_len = filename.length()+1;

      //Py_DebugFlag = 1;
      //Py_VerboseFlag = 1;

#ifdef PYTHON2
      char* sfilename = new char[fn_len];
      snprintf( sfilename, fn_len, "%s", filename.c_str() );
#else
      wchar_t* sfilename = new wchar_t[fn_len];
      swprintf( sfilename, fn_len, L"%s", filename.c_str() );
#endif
      if ( sfilename )
        Py_SetProgramName( sfilename );

      Py_InitializeEx( 1 );

      if ( sfilename )
        delete [] sfilename;
      if ( !Py_IsInitialized() )
        throw CG_FATAL( "PythonHandler" ) << "Failed to initialise the Python cards parser!";

      CG_INFO( "PythonHandler" )
        << "Initialised the Python cards parser\n\t"
        << "Python version: " << Py_GetVersion() << "\n\t"
        << "Platform: " << Py_GetPlatform() << ".";

      PyObject* cfg = PyImport_ImportModule( filename.c_str() ); // new
      if ( !cfg )
        throwPythonError( Form( "Failed to parse the configuration card %s", file ) );

      PyObject* process = PyObject_GetAttrString( cfg, PROCESS_NAME ); // new
      if ( !process )
        throwPythonError( Form( "Failed to extract a \"%s\" keyword from the configuration card %s", PROCESS_NAME, file ) );

      //--- list of process-specific parameters
      ParametersList proc_params;
      fillParameter( process, "processParameters", proc_params );

      //--- type of process to consider
      PyObject* pproc_name = getElement( process, MODULE_NAME ); // borrowed
      if ( !pproc_name )
        throwPythonError( Form( "Failed to extract the process name from the configuration card %s", file ) );
      const std::string proc_name = get<std::string>( pproc_name );

      //--- process mode
      params_.kinematics.mode = (KinematicsMode)proc_params.get<int>( "mode", (int)KinematicsMode::invalid );

      if ( proc_name == "lpair" )
        params_.setProcess( new Process::GamGamLL( proc_params ) );
      else if ( proc_name == "pptoll" || proc_name == "pptoff" )
        params_.setProcess( new Process::PPtoFF( proc_params ) );
      else if ( proc_name == "pptoww" )
        params_.setProcess( new Process::PPtoWW( proc_params ) );
      else {
        Process::generateFortranProcesses();
        for ( auto& proc : Process::FortranProcessesHandler::get().list() )
          if ( proc_name == std::string( proc.name ) )
            params_.setProcess( new Process::FortranKTProcess( proc_params, proc.name, proc.description, proc.method ) );
        if ( !params_.process() )
          throw CG_FATAL( "PythonHandler" ) << "Unrecognised process name: " << proc_name << "!";
      }

      //--- process kinematics
      PyObject* pin_kinematics = getElement( process, "inKinematics" ); // borrowed
      if ( pin_kinematics )
        parseIncomingKinematics( pin_kinematics );

      PyObject* pout_kinematics = getElement( process, "outKinematics" ); // borrowed
      if ( pout_kinematics )
        parseOutgoingKinematics( pout_kinematics );

      //--- taming functions
      PyObject* ptam = getElement( process, "tamingFunctions" ); // borrowed
      if ( ptam )
        parseTamingFunctions( ptam );

      Py_CLEAR( process );

      PyObject* plog = PyObject_GetAttrString( cfg, "logger" ); // new
      if ( plog ) {
        parseLogging( plog );
        Py_CLEAR( plog );
      }

      //--- hadroniser parameters
      PyObject* phad = PyObject_GetAttrString( cfg, "hadroniser" ); // new
      if ( phad ) {
        parseHadroniser( phad );
        Py_CLEAR( phad );
      }

      //--- generation parameters
      PyObject* pint = PyObject_GetAttrString( cfg, "integrator" ); // new
      if ( pint ) {
        parseIntegrator( pint );
        Py_CLEAR( pint );
      }

      PyObject* pgen = PyObject_GetAttrString( cfg, "generator" ); // new
      if ( pgen ) {
        parseGenerator( pgen );
        Py_CLEAR( pgen );
      }

      //--- finalisation
      Py_CLEAR( cfg );
    }

    PythonHandler::~PythonHandler()
    {
      if ( Py_IsInitialized() )
        Py_Finalize();
    }

    void
    PythonHandler::parseIncomingKinematics( PyObject* kin )
    {
      //--- retrieve the beams PDG ids
      std::vector<double> beams_pz;
      fillParameter( kin, "pz", beams_pz );
      if ( beams_pz.size() == 2 ) {
        params_.kinematics.incoming_beams.first.pz = beams_pz.at( 0 );
        params_.kinematics.incoming_beams.second.pz = beams_pz.at( 1 );
      }
      //--- retrieve the beams longitudinal momentum
      std::vector<int> beams_pdg;
      fillParameter( kin, "pdgIds", beams_pdg );
      if ( beams_pdg.size() == 2 ) {
        params_.kinematics.incoming_beams.first.pdg = (PDG)beams_pdg.at( 0 );
        params_.kinematics.incoming_beams.second.pdg = (PDG)beams_pdg.at( 1 );
      }
      double sqrt_s = -1.;
      fillParameter( kin, "cmEnergy", sqrt_s );
      fillParameter( kin, "kmrGridPath", params_.kinematics.kmr_grid_path );
      if ( sqrt_s != -1. )
        params_.kinematics.setSqrtS( sqrt_s );
      PyObject* psf = getElement( kin, "structureFunctions" ); // borrowed
      if ( psf )
        parseStructureFunctions( psf, params_.kinematics.structure_functions );
      std::vector<int> kt_fluxes;
      fillParameter( kin, "ktFluxes", kt_fluxes );
      if ( kt_fluxes.size() > 0 )
        params_.kinematics.incoming_beams.first.kt_flux = (KTFlux)kt_fluxes.at( 0 );
      if ( kt_fluxes.size() > 1 )
        params_.kinematics.incoming_beams.second.kt_flux = (KTFlux)kt_fluxes.at( 1 );
      std::vector<int> hi_beam1, hi_beam2;
      fillParameter( kin, "heavyIonA", hi_beam1 );
      if ( hi_beam1.size() == 2 )
        params_.kinematics.incoming_beams.first.pdg = HeavyIon{ (unsigned short)hi_beam1[0], (Element)hi_beam1[1] };
      fillParameter( kin, "heavyIonB", hi_beam2 );
      if ( hi_beam2.size() == 2 )
        params_.kinematics.incoming_beams.second.pdg = HeavyIon{ (unsigned short)hi_beam2[0], (Element)hi_beam2[1] };
    }

    void
    PythonHandler::parseStructureFunctions( PyObject* psf, std::shared_ptr<StructureFunctions>& sf_handler )
    {
      int str_fun = 0;
      fillParameter( psf, "id", str_fun );
      sf_handler = StructureFunctionsBuilder::get( (SF::Type)str_fun );
      switch( (SF::Type)str_fun ) {
        case SF::Type::LHAPDF: {
          auto sf = std::dynamic_pointer_cast<SF::LHAPDF>( params_.kinematics.structure_functions );
          fillParameter( psf, "pdfSet", sf->params.pdf_set );
          fillParameter( psf, "numFlavours", (unsigned int&)sf->params.num_flavours );
          fillParameter( psf, "pdfMember", (unsigned int&)sf->params.pdf_member );
          fillParameter( psf, "mode", (unsigned int&)sf->params.mode );
        } break;
        case SF::Type::MSTWgrid: {
          auto sf = std::dynamic_pointer_cast<mstw::Grid>( params_.kinematics.structure_functions );
          fillParameter( psf, "gridPath", sf->params.grid_path );
        } break;
        case SF::Type::Schaefer: {
          auto sf = std::dynamic_pointer_cast<SF::Schaefer>( params_.kinematics.structure_functions );
          fillParameter( psf, "Q2cut", sf->params.q2_cut );
          std::vector<double> w2_lims;
          fillParameter( psf, "W2limits", w2_lims );
          if ( w2_lims.size() != 0 ) {
            if ( w2_lims.size() != 2 )
              throwPythonError( Form( "Invalid size for W2limits attribute: %d != 2!", w2_lims.size() ) );
            else {
              sf->params.w2_lo = *std::min_element( w2_lims.begin(), w2_lims.end() );
              sf->params.w2_hi = *std::max_element( w2_lims.begin(), w2_lims.end() );
            }
          }
          PyObject* pcsf = getElement( psf, "continuumSF" ); // borrowed
          if ( pcsf )
            parseStructureFunctions( pcsf, sf->params.continuum_model );
          PyObject* ppsf = getElement( psf, "perturbativeSF" ); // borrowed
          if ( ppsf )
            parseStructureFunctions( ppsf, sf->params.perturbative_model );
          PyObject* prsf = getElement( psf, "resonancesSF" ); // borrowed
          if ( prsf )
            parseStructureFunctions( prsf, sf->params.resonances_model );
          fillParameter( psf, "higherTwist", (bool&)sf->params.higher_twist );
        } break;
        default: break;
      }
    }

    void
    PythonHandler::parseOutgoingKinematics( PyObject* kin )
    {
      PyObject* pparts = getElement( kin, "minFinalState" ); // borrowed
      if ( pparts && PyTuple_Check( pparts ) )
        for ( unsigned short i = 0; i < PyTuple_Size( pparts ); ++i )
          params_.kinematics.minimum_final_state.emplace_back( (PDG)get<int>( PyTuple_GetItem( pparts, i ) ) );

      PyObject* pcuts = getElement( kin, "cuts" ); // borrowed
      if ( pcuts )
        parseParticlesCuts( pcuts );

      // for LPAIR/collinear matrix elements
      fillLimits( kin, "q2", params_.kinematics.cuts.initial.q2 );

      // for the kT factorised matrix elements
      fillLimits( kin, "qt", params_.kinematics.cuts.initial.qt );
      fillLimits( kin, "phiqt", params_.kinematics.cuts.initial.phi_qt );
      fillLimits( kin, "ptdiff", params_.kinematics.cuts.central.pt_diff );
      fillLimits( kin, "phiptdiff", params_.kinematics.cuts.central.phi_pt_diff );
      fillLimits( kin, "rapiditydiff", params_.kinematics.cuts.central.rapidity_diff );

      // generic phase space limits
      fillLimits( kin, "rapidity", params_.kinematics.cuts.central.rapidity_single );
      fillLimits( kin, "eta", params_.kinematics.cuts.central.eta_single );
      fillLimits( kin, "pt", params_.kinematics.cuts.central.pt_single );

      fillLimits( kin, "ptsum", params_.kinematics.cuts.central.pt_sum );
      fillLimits( kin, "invmass", params_.kinematics.cuts.central.mass_sum );

      fillLimits( kin, "mx", params_.kinematics.cuts.remnants.mass_single );
    }

    void
    PythonHandler::parseParticlesCuts( PyObject* cuts )
    {
      if ( !PyDict_Check( cuts ) )
        throwPythonError( "Particle cuts object should be a dictionary!" );
      PyObject* pkey = nullptr, *pvalue = nullptr;
      Py_ssize_t pos = 0;
      while ( PyDict_Next( cuts, &pos, &pkey, &pvalue ) ) {
        const PDG pdg = (PDG)get<int>( pkey );
        fillLimits( pvalue, "pt", params_.kinematics.cuts.central_particles[pdg].pt_single );
        fillLimits( pvalue, "energy", params_.kinematics.cuts.central_particles[pdg].energy_single );
        fillLimits( pvalue, "eta", params_.kinematics.cuts.central_particles[pdg].eta_single );
        fillLimits( pvalue, "rapidity", params_.kinematics.cuts.central_particles[pdg].rapidity_single );
      }
    }

    void
    PythonHandler::parseLogging( PyObject* log )
    {
      fillParameter( log, "level", (int&)Logger::get().level );
      std::vector<std::string> enabled_modules;
      fillParameter( log, "enabledModules", enabled_modules );
      for ( const auto& mod : enabled_modules )
        Logger::get().addExceptionRule( mod );
    }

    void
    PythonHandler::parseIntegrator( PyObject* integr )
    {
      if ( !PyDict_Check( integr ) )
        throwPythonError( "Integrator object should be a dictionary!" );
      PyObject* palgo = getElement( integr, MODULE_NAME ); // borrowed
      if ( !palgo )
        throwPythonError( "Failed to retrieve the integration algorithm name!" );
      std::string algo = get<std::string>( palgo );
      if ( algo == "plain" )
        params_.integrator.type = Integrator::Type::plain;
      else if ( algo == "Vegas" ) {
        params_.integrator.type = Integrator::Type::Vegas;
        fillParameter( integr, "alpha", (double&)params_.integrator.vegas.alpha );
        fillParameter( integr, "iterations", params_.integrator.vegas.iterations );
        fillParameter( integr, "mode", (int&)params_.integrator.vegas.mode );
        fillParameter( integr, "verbosity", (int&)params_.integrator.vegas.verbose );
        std::string vegas_logging_output = "cerr";
        fillParameter( integr, "loggingOutput", vegas_logging_output );
        if ( vegas_logging_output == "cerr" )
          // redirect all debugging information to the error stream
          params_.integrator.vegas.ostream = stderr;
        else if ( vegas_logging_output == "cout" )
          // redirect all debugging information to the standard stream
          params_.integrator.vegas.ostream = stdout;
        else
          params_.integrator.vegas.ostream = fopen( vegas_logging_output.c_str(), "w" );
      }
      else if ( algo == "MISER" ) {
        params_.integrator.type = Integrator::Type::MISER;
        fillParameter( integr, "estimateFraction", (double&)params_.integrator.miser.estimate_frac );
        fillParameter( integr, "minCalls", params_.integrator.miser.min_calls );
        fillParameter( integr, "minCallsPerBisection", params_.integrator.miser.min_calls_per_bisection );
        fillParameter( integr, "alpha", (double&)params_.integrator.miser.alpha );
        fillParameter( integr, "dither", (double&)params_.integrator.miser.dither );
      }
      else
        throwPythonError( Form( "Invalid integration algorithm: %s", algo.c_str() ) );

      fillParameter( integr, "numFunctionCalls", params_.integrator.ncvg );
      fillParameter( integr, "seed", (unsigned long&)params_.integrator.rng_seed );
      unsigned int rng_engine;
      fillParameter( integr, "rngEngine", rng_engine );
      switch ( rng_engine ) {
        case 0: default: params_.integrator.rng_engine = (gsl_rng_type*)gsl_rng_mt19937; break;
        case 1: params_.integrator.rng_engine = (gsl_rng_type*)gsl_rng_taus2; break;
        case 2: params_.integrator.rng_engine = (gsl_rng_type*)gsl_rng_gfsr4; break;
        case 3: params_.integrator.rng_engine = (gsl_rng_type*)gsl_rng_ranlxs0; break;
      }
      fillParameter( integr, "chiSqCut", params_.integrator.vegas_chisq_cut );
    }

    void
    PythonHandler::parseGenerator( PyObject* gen )
    {
      if ( !PyDict_Check( gen ) )
        throwPythonError( "Generation information object should be a dictionary!" );
      params_.generation.enabled = true;
      fillParameter( gen, "treat", params_.generation.treat );
      fillParameter( gen, "numEvents", params_.generation.maxgen );
      fillParameter( gen, "printEvery", params_.generation.gen_print_every );
      fillParameter( gen, "numThreads", params_.generation.num_threads );
      fillParameter( gen, "numPoints", params_.generation.num_points );
    }

    void
    PythonHandler::parseTamingFunctions( PyObject* tf )
    {
      if ( !PyList_Check( tf ) )
        throwPythonError( "Taming functions list should be a list!" );

      for ( Py_ssize_t i = 0; i < PyList_Size( tf ); ++i ) {
        PyObject* pit = PyList_GetItem( tf, i ); // borrowed
        if ( !pit )
          continue;
        if ( !PyDict_Check( pit ) )
          throwPythonError( Form( "Item %d has invalid type %s", i, pit->ob_type->tp_name ) );
        PyObject* pvar = getElement( pit, "variable" ), *pexpr = getElement( pit, "expression" ); // borrowed
        params_.taming_functions->add( get<std::string>( pvar ).c_str(), get<std::string>( pexpr ).c_str() );
      }
    }

    void
    PythonHandler::parseHadroniser( PyObject* hadr )
    {
      if ( !PyDict_Check( hadr ) )
        throwPythonError( "Hadroniser object should be a dictionary!" );

      PyObject* pname = getElement( hadr, MODULE_NAME ); // borrowed
      if ( !pname )
        throwPythonError( "Hadroniser name is required!" );
      std::string hadr_name = get<std::string>( pname );

      fillParameter( hadr, "maxTrials", params_.hadroniser_max_trials );
      PyObject* pseed = getElement( hadr, "seed" ); // borrowed
      long long seed = -1ll;
      if ( pseed && is<int>( pseed ) ) {
        seed = PyLong_AsLongLong( pseed );
        CG_DEBUG( "PythonHandler:hadroniser" ) << "Hadroniser seed set to " << seed;
      }
      if ( hadr_name == "pythia8" ) {
        params_.setHadroniser( new Hadroniser::Pythia8Hadroniser( params_ ) );
        std::vector<std::string> config;
        auto pythia8 = dynamic_cast<Hadroniser::Pythia8Hadroniser*>( params_.hadroniser() );
        pythia8->setSeed( seed );
        fillParameter( hadr, "pythiaPreConfiguration", config );
        pythia8->readStrings( config );
        pythia8->init();
        fillParameter( hadr, "pythiaConfiguration", config );
        pythia8->readStrings( config );
        fillParameter( hadr, "pythiaProcessConfiguration", config );
        pythia8->readStrings( config );
      }
    }
  }
}

#endif

