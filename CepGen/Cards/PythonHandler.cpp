#include "PythonHandler.h"
#include "CepGen/Core/Exception.h"

#ifdef PYTHON

#include "CepGen/Core/TamingFunction.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoFF.h"
#include "CepGen/Processes/PPtoWW.h"
#include "CepGen/Processes/FortranKTProcess.h"

#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/StructureFunctions/GenericLHAPDF.h"
#include "CepGen/StructureFunctions/MSTWGrid.h"
#include "CepGen/StructureFunctions/Schaefer.h"

#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

#include <algorithm>
#include <frameobject.h>

extern "C"
{
  extern void nucl_to_ff_( double& );
}

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

      //--- type of process to consider
      PyObject* pproc_name = getElement( process, MODULE_NAME ); // borrowed
      if ( !pproc_name )
        throwPythonError( Form( "Failed to extract the process name from the configuration card %s", file ) );
      const std::string proc_name = decode( pproc_name );

      if ( proc_name == "lpair" )
        params_.setProcess( new Process::GamGamLL );
      else if ( proc_name == "pptoll" || proc_name == "pptoff" )
        params_.setProcess( new Process::PPtoFF );
      else if ( proc_name == "pptoww" )
        params_.setProcess( new Process::PPtoWW );
      else if ( proc_name == "patoll" )
        params_.setProcess( new Process::FortranKTProcess( "nucltoff", "(p/A)(p/A) ↝ (g/ɣ)ɣ → f⁺f¯", nucl_to_ff_ ) );
      else throw CG_FATAL( "PythonHandler" ) << "Unrecognised process: " << proc_name << ".";

      //--- process mode
      fillParameter( process, "mode", (int&)params_.kinematics.mode );

      //--- process kinematics
      PyObject* pin_kinematics = getElement( process, "inKinematics" ); // borrowed
      if ( pin_kinematics )
        parseIncomingKinematics( pin_kinematics );

      PyObject* pout_kinematics = getElement( process, "outKinematics" ); // borrowed
      if ( pout_kinematics )
        parseOutgoingKinematics( pout_kinematics );

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

      //--- taming functions
      PyObject* ptam = PyObject_GetAttrString( cfg, "tamingFunctions" ); // new
      if ( ptam ) {
        parseTamingFunctions( ptam );
        Py_CLEAR( ptam );
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
      PyObject* ppz = getElement( kin, "pz" ); // borrowed
      if ( ppz && PyTuple_Check( ppz ) && PyTuple_Size( ppz ) == 2 ) {
        double pz0 = PyFloat_AsDouble( PyTuple_GetItem( ppz, 0 ) );
        double pz1 = PyFloat_AsDouble( PyTuple_GetItem( ppz, 1 ) );
        params_.kinematics.incoming_beams.first.pz = pz0;
        params_.kinematics.incoming_beams.second.pz = pz1;
      }
      double sqrt_s = -1.;
      fillParameter( kin, "cmEnergy", sqrt_s );
      fillParameter( kin, "kmrGridPath", params_.kinematics.kmr_grid_path );
      if ( sqrt_s != -1. )
        params_.kinematics.setSqrtS( sqrt_s );
      PyObject* psf = getElement( kin, "structureFunctions" ); // borrowed
      if ( psf )
        parseStructureFunctions( psf );
      std::vector<int> kt_fluxes;
      fillParameter( kin, "ktFluxes", kt_fluxes );
      if ( kt_fluxes.size() > 0 )
        params_.kinematics.incoming_beams.first.kt_flux = kt_fluxes.at( 0 );
      if ( kt_fluxes.size() > 1 )
        params_.kinematics.incoming_beams.second.kt_flux = kt_fluxes.at( 1 );
      std::vector<int> hi_beam1, hi_beam2;
      fillParameter( kin, "heavyIonA", hi_beam1 );
      if ( hi_beam1.size() == 2 )
        params_.kinematics.incoming_beams.first.hi = Kinematics::HeavyIon{ (unsigned short)hi_beam1[0], (unsigned short)hi_beam1[1] };
      fillParameter( kin, "heavyIonB", hi_beam2 );
      if ( hi_beam2.size() == 2 )
        params_.kinematics.incoming_beams.second.hi = Kinematics::HeavyIon{ (unsigned short)hi_beam2[0], (unsigned short)hi_beam2[1] };
    }

    void
    PythonHandler::parseStructureFunctions( PyObject* psf )
    {
      int str_fun = 0;
      fillParameter( psf, "id", str_fun );
      params_.kinematics.structure_functions.reset( StructureFunctionsBuilder::get( (SF::Type)str_fun ) );
      switch( (SF::Type)str_fun ) {
        case SF::Type::GenericLHAPDF: {
          auto sf = dynamic_cast<SF::GenericLHAPDF*>( params_.kinematics.structure_functions.get() );
          fillParameter( psf, "pdfSet", sf->params.pdf_set );
          fillParameter( psf, "numFlavours", (unsigned int&)sf->params.num_flavours );
        } break;
        case SF::Type::MSTWgrid: {
          auto sf = dynamic_cast<MSTW::Grid*>( params_.kinematics.structure_functions.get() );
          fillParameter( psf, "gridPath", sf->params.grid_path );
        } break;
        case SF::Type::Schaefer: {
          auto sf = dynamic_cast<SF::Schaefer*>( params_.kinematics.structure_functions.get() );
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
            fillParameter( pcsf, "id", sf->params.cont_model );
          PyObject* prsf = getElement( psf, "resonancesSF" ); // borrowed
          if ( prsf )
            fillParameter( prsf, "id", sf->params.res_model );
          fillParameter( psf, "higherTwist", (bool&)sf->params.higher_twist );
        } break;
        default: break;
      }
    }

    void
    PythonHandler::parseOutgoingKinematics( PyObject* kin )
    {
      PyObject* ppair = getElement( kin, "pair" ); // borrowed
      if ( ppair ) {
        if ( isInteger( ppair ) ) {
          PDG pair = (PDG)asInteger( ppair );
          params_.kinematics.central_system = { pair, pair };
        }
        else if ( PyTuple_Check( ppair ) ) {
          if ( PyTuple_Size( ppair ) != 2 )
            throw CG_FATAL( "PythonHandler" ) << "Invalid value for in_kinematics.pair!";
          PDG pair1 = (PDG)asInteger( PyTuple_GetItem( ppair, 0 ) );
          PDG pair2 = (PDG)asInteger( PyTuple_GetItem( ppair, 1 ) );
          params_.kinematics.central_system = { pair1, pair2 };
        }
      }

      PyObject* pparts = getElement( kin, "minFinalState" ); // borrowed
      if ( pparts && PyTuple_Check( pparts ) )
        for ( unsigned short i = 0; i < PyTuple_Size( pparts ); ++i )
          params_.kinematics.minimum_final_state.emplace_back( (PDG)asInteger( PyTuple_GetItem( pparts, i ) ) );

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
        const PDG pdg = (PDG)asInteger( pkey );
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
      std::string algo = decode( palgo );
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
        params_.taming_functions->add( decode( pvar ).c_str(), decode( pexpr ).c_str() );
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
      std::string hadr_name = decode( pname );

      fillParameter( hadr, "maxTrials", params_.hadroniser_max_trials );
      PyObject* pseed = getElement( hadr, "seed" ); // borrowed
      long long seed = -1ll;
      if ( pseed && isInteger( pseed ) ) {
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

    //------------------------------------------------------------------
    // Python API helpers
    //------------------------------------------------------------------

    std::string
    PythonHandler::getPythonPath( const char* file )
    {
      std::string s_filename = file;
      s_filename = s_filename.substr( 0, s_filename.find_last_of( "." ) ); // remove the extension
      std::replace( s_filename.begin(), s_filename.end(), '/', '.' ); // replace all '/' by '.'
      return s_filename;
    }

    void
    PythonHandler::throwPythonError( const std::string& message )
    {
      PyObject* ptype = nullptr, *pvalue = nullptr, *ptraceback_obj = nullptr;
      // retrieve error indicator and clear it to handle ourself the error
      PyErr_Fetch( &ptype, &pvalue, &ptraceback_obj );
      PyErr_Clear();
      // ensure the objects retrieved are properly normalised and point to compatible objects
      PyErr_NormalizeException( &ptype, &pvalue, &ptraceback_obj );
      std::ostringstream oss; oss << message;
      if ( ptype != nullptr ) { // we can start the traceback
        oss << "\n\tError: "
#ifdef PYTHON2
            << PyString_AsString( PyObject_Str( pvalue ) ); // deprecated in python v3+
#else
            << PyUnicode_AsUTF8( PyObject_Str( pvalue ) );
#endif
        PyTracebackObject* ptraceback = (PyTracebackObject*)ptraceback_obj;
        string tabul = "↪ ";
        if ( ptraceback != nullptr ) {
          while ( ptraceback->tb_next != nullptr ) {
            PyFrameObject* pframe = ptraceback->tb_frame;
            if ( pframe != nullptr ) {
              int line = PyCode_Addr2Line( pframe->f_code, pframe->f_lasti );
#ifdef PYTHON2
              const char* filename = PyString_AsString( pframe->f_code->co_filename );
              const char* funcname = PyString_AsString( pframe->f_code->co_name );
#else
              const char* filename = PyUnicode_AsUTF8( pframe->f_code->co_filename );
              const char* funcname = PyUnicode_AsUTF8( pframe->f_code->co_name );
#endif
              oss << Form( "\n\t%s%s on %s (line %d)", tabul.c_str(), boldify( funcname ).c_str(), filename, line );
            }
            else
              oss << Form( "\n\t%s issue in line %d", tabul.c_str(), ptraceback->tb_lineno );
            tabul = string( "  " )+tabul;
            ptraceback = ptraceback->tb_next;
          }
        }
      }
      Py_Finalize();
      throw CG_FATAL( "PythonHandler:error" ) << oss.str();
    }

    std::string
    PythonHandler::decode( PyObject* obj )
    {
      std::string out;
#ifdef PYTHON2
      out = PyString_AsString( obj ); // deprecated in python v3+
#else
      PyObject* pstr = PyUnicode_AsEncodedString( obj, "utf-8", "strict" ); // new
      if ( !pstr )
        throwPythonError( "Failed to decode a Python object!" );
      out = PyBytes_AS_STRING( pstr );
      Py_CLEAR( pstr );
#endif
      return out;
    }

    PyObject*
    PythonHandler::encode( const char* str )
    {
      PyObject* obj = PyUnicode_FromString( str ); // new
      if ( !obj )
        throwPythonError( Form( "Failed to encode the following string:\n\t%s", str ) );
      return obj;
    }

    PyObject*
    PythonHandler::getElement( PyObject* obj, const char* key )
    {
      PyObject* pout = nullptr, *nink = encode( key );
      if ( !nink )
        return pout;
      pout = PyDict_GetItem( obj, nink ); // borrowed
      Py_CLEAR( nink );
      if ( pout )
        CG_DEBUG( "PythonHandler:getElement" )
          << "retrieved " << pout->ob_type->tp_name << " element \"" << key << "\" "
          << "from " << obj->ob_type->tp_name << " object\n\t"
          << "new reference count: " << pout->ob_refcnt;
      else
        CG_DEBUG( "PythonHandler:getElement" )
          << "did not retrieve a valid element \"" << key << "\"";
      return pout;
    }

    void
    PythonHandler::fillLimits( PyObject* obj, const char* key, Limits& lim )
    {
      PyObject* pobj = getElement( obj, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) )
        throw CG_FATAL( "PythonHandler:fillLimits" ) << "Invalid value retrieved for " << key << ".";
      if ( PyTuple_Size( pobj ) < 1 )
        throw CG_FATAL( "PythonHandler:fillLimits" ) << "Invalid number of values unpacked for " << key << "!";
      double min = PyFloat_AsDouble( PyTuple_GetItem( pobj, 0 ) );
      lim.min() = min;
      if ( PyTuple_Size( pobj ) > 1 ) {
        double max = PyFloat_AsDouble( PyTuple_GetItem( pobj, 1 ) );
        if ( max != -1 )
          lim.max() = max;
      }
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, bool& out )
    {
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyBool_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      fillParameter( parent, key, (int&)out );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, int& out )
    {
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
#ifdef PYTHON2
      if ( !PyInt_Check( pobj ) && !PyBool_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      out = PyInt_AsLong( pobj );
#else
      if ( !PyLong_Check( pobj ) && !PyBool_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      out = PyLong_AsLong( pobj );
#endif
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, unsigned long& out )
    {
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyLong_Check( pobj )
#ifdef PYTHON2
        && !PyInt_Check( pobj )
#endif
      )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      if ( PyLong_Check( pobj ) )
        out = PyLong_AsUnsignedLong( pobj );
#ifdef PYTHON2
      else if ( PyInt_Check( pobj ) )
        out = PyInt_AsUnsignedLongMask( pobj );
#endif
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, unsigned int& out )
    {
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
#ifdef PYTHON2
      if ( !PyInt_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      out = PyInt_AsUnsignedLongMask( pobj );
#else
      if ( !PyLong_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      out = PyLong_AsUnsignedLong( pobj );
#endif
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, double& out )
    {
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyFloat_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      out = PyFloat_AsDouble( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, std::string& out )
    {
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !
#ifdef PYTHON2
        PyString_Check( pobj )
#else
        PyUnicode_Check( pobj )
#endif
      )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      out = decode( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, std::vector<double>& out )
    {
      out.clear();
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      for ( Py_ssize_t i = 0; i < PyTuple_Size( pobj ); ++i ) {
        PyObject* pit = PyTuple_GetItem( pobj, i ); // borrowed
        if ( !PyFloat_Check( pit ) )
          continue;
        out.emplace_back( PyFloat_AsDouble( pit ) );
      }
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, std::vector<std::string>& out )
    {
      out.clear();
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      for ( Py_ssize_t i = 0; i < PyTuple_Size( pobj ); ++i ) {
        PyObject* pit = PyTuple_GetItem( pobj, i ); // borrowed
        out.emplace_back( decode( pit ) );
      }
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, std::vector<int>& out )
    {
      out.clear();
      PyObject* pobj = getElement( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      for ( Py_ssize_t i = 0; i < PyTuple_Size( pobj ); ++i ) {
        PyObject* pit = PyTuple_GetItem( pobj, i );
        if ( !isInteger( pit ) )
          throwPythonError( Form( "Object %d has invalid type", i ) );
        out.emplace_back( asInteger( pit ) );
      }
    }

    bool
    PythonHandler::isInteger( PyObject* obj )
    {
#ifdef PYTHON2
      return PyInt_Check( obj );
#else
      return PyLong_Check( obj );
#endif
    }

    int
    PythonHandler::asInteger( PyObject* obj )
    {
#ifdef PYTHON2
      return PyInt_AsLong( obj );
#else
      return PyLong_AsLong( obj );
#endif
    }
  }
}

#endif

