#include "PythonHandler.h"
#include "CepGen/Core/Exception.h"

#ifdef PYTHON

#include "CepGen/Core/TamingFunction.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"
#include "CepGen/Processes/PPtoWW.h"

#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

#include <algorithm>
#include <frameobject.h>

#ifdef PYTHIA8
void
feedPythia( CepGen::Hadroniser::Pythia8Hadroniser* py8, PyObject* hadr, const char* config )
{
  PyObject* ppc = CepGen::Cards::PythonHandler::getElement( hadr, config );
  if ( !ppc )
    return;
  if ( !PyTuple_Check( ppc ) ) {
    Py_DECREF( ppc );
    return;
  }
  for ( Py_ssize_t i = 0; i < PyTuple_Size( ppc ); ++i ) {
    PyObject* pln = PyTuple_GetItem( ppc, i );
    if ( !pln )
      continue;
    py8->readString( CepGen::Cards::PythonHandler::decode( pln ) );
  }
}
#endif

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for CepGen input cards
    PythonHandler::PythonHandler( const char* file ) :
      sfilename_( nullptr ), cfg_( nullptr )
    {
      setenv( "PYTHONPATH", ".:..:Cards", 1 );
      std::string filename = getPythonPath( file );
      const size_t fn_len = filename.length()+1;
#ifdef PYTHON2
      sfilename_ = new char[fn_len];
      snprintf( sfilename_, fn_len, "%s", filename.c_str() );
#else
      sfilename_ = new wchar_t[fn_len];
      swprintf( sfilename_, fn_len, L"%s", filename.c_str() );
#endif
      Py_SetProgramName( sfilename_ );

      Py_InitializeEx( 0 );
      if ( !Py_IsInitialized() )
        throw CG_FATAL( "PythonHandler" ) << "Failed to initialise the python parser!";

      CG_DEBUG( "PythonHandler" )
        << "Initialised the Python cards parser\n\t"
        << "Python version: " << Py_GetVersion() << "\n\t"
        << "Platform: " << Py_GetPlatform() << ".";

      PyObject* fn = encode( filename.c_str() );
      if ( !fn )
        throwPythonError( Form( "Failed to encode the configuration filename %s", filename.c_str() ) );

      cfg_ = PyImport_Import( fn );
      Py_DECREF( fn );
      if ( !cfg_ )
        throwPythonError( Form( "Failed to parse the configuration card %s", file ) );

      PyObject* process = PyObject_GetAttrString( cfg_, "process" );
      if ( !process )
        throwPythonError( Form( "Failed to extract a \"process\" keyword from the configuration card %s", file ) );

      //--- type of process to consider
      PyObject* pproc_name = PyDict_GetItem( process, encode( module_name_ ) );
      if ( !pproc_name )
        throwPythonError( Form( "Failed to extract the process name from the configuration card %s", file ) );

      const std::string proc_name = decode( pproc_name );
      if ( proc_name == "lpair" )
        params_.setProcess( new Process::GamGamLL );
      else if ( proc_name == "pptoll" )
        params_.setProcess( new Process::PPtoLL );
      else if ( proc_name == "pptoww" )
        params_.setProcess( new Process::PPtoWW );
      else throw CG_FATAL( "PythonHandler" ) << "Unrecognised process: " << proc_name << ".";

      //--- process mode
      getParameter( process, "mode", (int&)params_.kinematics.mode );

      //--- process kinematics
      PyObject* pin_kinematics = getElement( process, "inKinematics" );
      if ( pin_kinematics )
        parseIncomingKinematics( pin_kinematics );

      PyObject* pout_kinematics = getElement( process, "outKinematics" );
      if ( pout_kinematics )
        parseOutgoingKinematics( pout_kinematics );

      Py_DECREF( process );

      PyObject* plog = PyObject_GetAttrString( cfg_, "logger" );
      if ( plog ) {
        parseLogging( plog );
        Py_DECREF( plog );
      }

      //--- hadroniser parameters
      PyObject* phad = PyObject_GetAttrString( cfg_, "hadroniser" );
      if ( phad ) {
        parseHadroniser( phad );
        Py_DECREF( phad );
      }

      //--- generation parameters
      PyObject* pint = PyObject_GetAttrString( cfg_, "integrator" );
      if ( pint ) {
        parseIntegrator( pint );
        Py_DECREF( pint );
      }

      PyObject* pgen = PyObject_GetAttrString( cfg_, "generator" );
      if ( pgen ) {
        parseGenerator( pgen );
        Py_DECREF( pgen );
      }

      //--- taming functions
      PyObject* ptam = PyObject_GetAttrString( cfg_, "tamingFunctions" );
      if ( ptam ) {
        parseTamingFunctions( ptam );
        Py_DECREF( ptam );
      }
    }

    PythonHandler::~PythonHandler()
    {
      Py_DECREF( cfg_ );
#ifdef PYTHON2
      Py_Finalize();
#else
      if ( Py_IsInitialized() || Py_FinalizeEx() != 0 )
        throw CG_FATAL( "PythonHandler" ) << "Failed to unregister the python parser!";
#endif

      if ( sfilename_ )
        delete [] sfilename_;
    }

    void
    PythonHandler::parseIncomingKinematics( PyObject* kin )
    {
      PyObject* ppz = getElement( kin, "pz" );
      if ( ppz ) {
        if ( PyTuple_Check( ppz ) && PyTuple_Size( ppz ) == 2 ) {
          double pz0 = PyFloat_AsDouble( PyTuple_GetItem( ppz, 0 ) );
          double pz1 = PyFloat_AsDouble( PyTuple_GetItem( ppz, 1 ) );
          params_.kinematics.inp = { pz0, pz1 };
        }
        Py_DECREF( ppz );
      }
      double sqrt_s = -1.;
      getParameter( kin, "cmEnergy", sqrt_s );
      if ( sqrt_s != -1. )
        params_.kinematics.setSqrtS( sqrt_s );
      getParameter( kin, "structureFunctions", (int&)params_.kinematics.structure_functions );
    }

    void
    PythonHandler::parseOutgoingKinematics( PyObject* kin )
    {
      PyObject* ppair = getElement( kin, "pair" );
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
        Py_DECREF( ppair );
      }

      PyObject* pparts = getElement( kin, "minFinalState" );
      if ( pparts ) {
        if ( PyTuple_Check( pparts ) )
          for ( unsigned short i = 0; i < PyTuple_Size( pparts ); ++i )
            params_.kinematics.minimum_final_state.emplace_back( (PDG)asInteger( PyTuple_GetItem( pparts, i ) ) );
        Py_DECREF( pparts );
      }

      PyObject* pcuts = getElement( kin, "cuts" );
      if ( pcuts && PyDict_Check( pcuts ) ) parseParticlesCuts( pcuts );

      // for LPAIR/collinear matrix elements
      getLimits( kin, "q2", params_.kinematics.cuts.initial[Cuts::q2] );

      // for the kT factorised matrix elements
      getLimits( kin, "qt", params_.kinematics.cuts.initial[Cuts::qt] );
      getLimits( kin, "phiqt", params_.kinematics.cuts.initial[Cuts::phi_qt] );
      getLimits( kin, "ptdiff", params_.kinematics.cuts.central[Cuts::pt_diff] );
      getLimits( kin, "phiptdiff", params_.kinematics.cuts.central[Cuts::phi_pt_diff] );
      getLimits( kin, "rapiditydiff", params_.kinematics.cuts.central[Cuts::rapidity_diff] );

      // generic phase space limits
      getLimits( kin, "rapidity", params_.kinematics.cuts.central[Cuts::rapidity_single] );
      getLimits( kin, "eta", params_.kinematics.cuts.central[Cuts::eta_single] );
      getLimits( kin, "pt", params_.kinematics.cuts.central[Cuts::pt_single] );

      getLimits( kin, "mx", params_.kinematics.cuts.remnants[Cuts::mass_single] );
    }

    void
    PythonHandler::parseParticlesCuts( PyObject* cuts )
    {
      PyObject* pkey = nullptr, *pvalue = nullptr;
      Py_ssize_t pos = 0;
      while ( PyDict_Next( cuts, &pos, &pkey, &pvalue ) ) {
        PDG pdg = (PDG)asInteger( pkey );
        getLimits( pvalue, "pt", params_.kinematics.cuts.central_particles[pdg][Cuts::pt_single] );
        getLimits( pvalue, "energy", params_.kinematics.cuts.central_particles[pdg][Cuts::energy_single] );
        getLimits( pvalue, "eta", params_.kinematics.cuts.central_particles[pdg][Cuts::eta_single] );
        getLimits( pvalue, "rapidity", params_.kinematics.cuts.central_particles[pdg][Cuts::rapidity_single] );
      }
    }

    void
    PythonHandler::parseLogging( PyObject* log )
    {
      getParameter( log, "level", (int&)Logger::get().level );
      std::vector<std::string> enabled_modules;
      getParameter( log, "enabledModules", enabled_modules );
      for ( const auto& mod : enabled_modules )
        Logger::get().addExceptionRule( mod );
    }

    void
    PythonHandler::parseIntegrator( PyObject* integr )
    {
      if ( !PyDict_Check( integr ) )
        throwPythonError( "Integrator object should be a dictionary!" );
      PyObject* palgo = getElement( integr, module_name_ );
      if ( !palgo )
        throwPythonError( "Failed to retrieve the integration algorithm name!" );
      std::string algo = decode( palgo );
      if ( algo == "plain" )
        params_.integrator.type = Integrator::Type::plain;
      else if ( algo == "Vegas" ) {
        params_.integrator.type = Integrator::Type::Vegas;
        getParameter( integr, "alpha", (double&)params_.integrator.vegas.alpha );
        getParameter( integr, "iterations", params_.integrator.vegas.iterations );
        getParameter( integr, "mode", (int&)params_.integrator.vegas.mode );
        getParameter( integr, "verbosity", (int&)params_.integrator.vegas.verbose );
        std::string vegas_logging_output = "cerr";
        getParameter( integr, "loggingOutput", vegas_logging_output );
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
        getParameter( integr, "estimateFraction", (double&)params_.integrator.miser.estimate_frac );
        getParameter( integr, "minCalls", params_.integrator.miser.min_calls );
        getParameter( integr, "minCallsPerBisection", params_.integrator.miser.min_calls_per_bisection );
        getParameter( integr, "alpha", (double&)params_.integrator.miser.alpha );
        getParameter( integr, "dither", (double&)params_.integrator.miser.dither );
      }
      else
        throwPythonError( Form( "Invalid integration algorithm: %s", algo.c_str() ) );

      getParameter( integr, "numPoints", params_.integrator.npoints );
      getParameter( integr, "numFunctionCalls", params_.integrator.ncvg );
      getParameter( integr, "seed", (unsigned long&)params_.integrator.rng_seed );
      unsigned int rng_engine;
      getParameter( integr, "rngEngine", rng_engine );
      switch ( rng_engine ) {
        case 0: default: params_.integrator.rng_engine = (gsl_rng_type*)gsl_rng_mt19937; break;
        case 1: params_.integrator.rng_engine = (gsl_rng_type*)gsl_rng_taus2; break;
        case 2: params_.integrator.rng_engine = (gsl_rng_type*)gsl_rng_gfsr4; break;
        case 3: params_.integrator.rng_engine = (gsl_rng_type*)gsl_rng_ranlxs0; break;
      }
      getParameter( integr, "chiSqCut", params_.integrator.vegas_chisq_cut );
    }

    void
    PythonHandler::parseGenerator( PyObject* gen )
    {
      if ( !PyDict_Check( gen ) )
        throwPythonError( "Generation information object should be a dictionary!" );
      params_.generation.enabled = true;
      getParameter( gen, "numEvents", params_.generation.maxgen );
      getParameter( gen, "printEvery", params_.generation.gen_print_every );
      getParameter( gen, "numThreads", params_.generation.num_threads );
    }

    void
    PythonHandler::parseTamingFunctions( PyObject* tf )
    {
      if ( !PyList_Check( tf ) )
        throwPythonError( "Taming functions list should be a list!" );

      for ( Py_ssize_t i = 0; i < PyList_Size( tf ); ++i ) {
        PyObject* pit = PyList_GetItem( tf, i );
        if ( !PyDict_Check( pit ) )
          throwPythonError( Form( "Item %d is invalid", i ) );
        PyObject* pvar = getElement( pit, "variable" );
        PyObject* pexpr = getElement( pit, "expression" );
        params_.taming_functions->add( decode( pvar ).c_str(), decode( pexpr ).c_str() );
        Py_DECREF( pvar );
        Py_DECREF( pexpr );
      }
    }

    void
    PythonHandler::parseHadroniser( PyObject* hadr )
    {
      if ( !PyDict_Check( hadr ) )
        throwPythonError( "Hadroniser object should be a dictionary!" );

      PyObject* pname = getElement( hadr, module_name_ );
      if ( !pname )
        throwPythonError( "Hadroniser name is required!" );

      std::string hadr_name = decode( pname );

      if ( hadr_name == "pythia8" ) {
        getParameter( hadr, "maxTrials", params_.hadroniser_max_trials );
#ifdef PYTHIA8
        params_.setHadroniser( new Hadroniser::Pythia8Hadroniser( params_ ) );
        auto pythia8 = dynamic_cast<Hadroniser::Pythia8Hadroniser*>( params_.hadroniser() );
        PyObject* pseed = getElement( hadr, "seed" );
        long long seed = -1ll;
        if ( pseed ) {
          if ( isInteger( pseed ) )
            seed = PyLong_AsLongLong( pseed );
          Py_DECREF( pseed );
        }
        pythia8->setSeed( seed );
        feedPythia( pythia8, hadr, "pythiaPreConfiguration" );
        pythia8->init();
        feedPythia( pythia8, hadr, "pythiaConfiguration" );
        feedPythia( pythia8, hadr, "pythiaProcessConfiguration" );
#else
        CG_Warning( "PythonHandler" )
          << "Pythia8 is not linked to this instance... "
          << "Ignoring this part of the configuration file.";
#endif
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
      //size_t python_subdir_pos = s_filename.find( "/python" );
      //if ( s_filename.size() > 7 && python_subdir_pos != std::string::npos )
      //  s_filename.erase( python_subdir_pos, python_subdir_pos+7 );
      std::replace( s_filename.begin(), s_filename.end(), '/', '.' ); // replace all '/' by '.'
      return s_filename;
    }

    void
    PythonHandler::throwPythonError( const std::string& message )
    {
      PyObject* ptype = nullptr, *pvalue = nullptr, *ptraceback_obj = nullptr;
      PyErr_Fetch( &ptype, &pvalue, &ptraceback_obj );
      PyErr_Clear();
      PyErr_NormalizeException( &ptype, &pvalue, &ptraceback_obj );
      std::ostringstream oss; oss << message;
      if ( ptype == nullptr ) {
        Py_Finalize();
        throw CG_FATAL( "PythonHandler" ) << oss.str();
      }

      oss << "\n\tError: "
#ifdef PYTHON2
          << PyString_AsString( PyObject_Str( pvalue ) ); // deprecated in python v3+
#else
          << _PyUnicode_AsString( PyObject_Str( pvalue ) );
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
            const char* filename = _PyUnicode_AsString( pframe->f_code->co_filename );
            const char* funcname = _PyUnicode_AsString( pframe->f_code->co_name );
#endif
            oss << Form( "\n\t%s%s on %s (line %d)", tabul.c_str(), boldify( funcname ).c_str(), filename, line );
            tabul = string( "  " )+tabul;
          }
          else
            oss << Form( "\n\t\tissue in line %d", ptraceback->tb_lineno );
          ptraceback = ptraceback->tb_next;
        }
      }
      /*PyThreadState* ptstate = PyThreadState_GET();
      if ( ptstate != nullptr && ptstate->frame != nullptr ) {
        PyFrameObject* pframe = ptstate->frame;
        while ( pframe != nullptr ) {
          int line = PyCode_Addr2Line( pframe->f_code, pframe->f_lasti );
#ifdef PYTHON2
          const char* filename = PyString_AsString( pframe->f_code->co_filename );
          const char* funcname = PyString_AsString( pframe->f_code->co_name );
#else
          const char* filename = _PyUnicode_AsString( pframe->f_code->co_filename );
          const char* funcname = _PyUnicode_AsString( pframe->f_code->co_name );
#endif
          oss << Form( "Stack trace:\n\t\t%s(%d): %s\n", filename, line, funcname );
          pframe = pframe->f_back;
        }
      }*/
      Py_Finalize();
      throw CG_FATAL( "PythonError" ) << oss.str();
    }

    std::string
    PythonHandler::decode( PyObject* obj )
    {
#ifdef PYTHON2
      const char* str = PyString_AsString( obj ); // deprecated in python v3+
#else
      const char* str = _PyUnicode_AsString( obj );
#endif
      //Py_DECREF( obj );
      if ( !str )
        throwPythonError( "Failed to decode a Python object!" );
      return std::string( str );
    }

    PyObject*
    PythonHandler::encode( const char* str )
    {
      PyObject* obj = PyUnicode_FromString( str );
      if ( !obj )
        throwPythonError( Form( "Failed to encode the following string:\n\t%s", str ) );
      return obj;
    }

    PyObject*
    PythonHandler::getElement( PyObject* obj, const char* key )
    {
      PyObject* pout = nullptr;
      PyObject* nink = encode( key );
      if ( !nink )
        return pout;
      pout = PyDict_GetItem( obj, nink );
      Py_DECREF( nink );
      return pout;
    }

    void
    PythonHandler::getLimits( PyObject* obj, const char* key, Limits& lim )
    {
      PyObject* pobj = getElement( obj, key );
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) )
        throw CG_FATAL( "PythonHandler" ) << "Invalid value retrieved for " << key << ".";
      if ( PyTuple_Size( pobj ) < 1 )
        throw CG_FATAL( "PythonHandler" ) << "Invalid number of values unpacked for " << key << "!";
      double min = PyFloat_AsDouble( PyTuple_GetItem( pobj, 0 ) );
      lim.min() = min;
      if ( PyTuple_Size( pobj ) > 1 ) {
        double max = PyFloat_AsDouble( PyTuple_GetItem( pobj, 1 ) );
        if ( max != -1 )
          lim.max() = max;
      }
      Py_DECREF( pobj );
    }

    void
    PythonHandler::getParameter( PyObject* parent, const char* key, int& out )
    {
      PyObject* pobj = getElement( parent, key );
      if ( !pobj )
        return;
#ifdef PYTHON2
      if ( !PyInt_Check( pobj ) ) {
        Py_DECREF( pobj );
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      }
      out = PyInt_AsLong( pobj );
#else
      if ( !PyLong_Check( pobj ) ) {
        Py_DECREF( pobj );
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      }
      out = PyLong_AsLong( pobj );
#endif
      Py_DECREF( pobj );
    }

    void
    PythonHandler::getParameter( PyObject* parent, const char* key, unsigned long& out )
    {
      PyObject* pobj = getElement( parent, key );
      if ( !pobj )
        return;
      if ( !PyLong_Check( pobj )
#ifdef PYTHON2
        && !PyInt_Check( pobj )
#endif
      ) {
        Py_DECREF( pobj );
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      }
      if ( PyLong_Check( pobj ) )
        out = PyLong_AsUnsignedLong( pobj );
#ifdef PYTHON2
      else if ( PyInt_Check( pobj ) )
        out = PyInt_AsUnsignedLongMask( pobj );
#endif
      Py_DECREF( pobj );
    }

    void
    PythonHandler::getParameter( PyObject* parent, const char* key, unsigned int& out )
    {
      PyObject* pobj = getElement( parent, key );
      if ( !pobj )
        return;
#ifdef PYTHON2
      if ( !PyInt_Check( pobj ) ) {
        Py_DECREF( pobj );
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      }
      out = PyInt_AsUnsignedLongMask( pobj );
#else
      if ( !PyLong_Check( pobj ) ) {
        Py_DECREF( pobj );
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      }
      out = PyLong_AsUnsignedLong( pobj );
#endif
      Py_DECREF( pobj );
    }

    void
    PythonHandler::getParameter( PyObject* parent, const char* key, double& out )
    {
      PyObject* pobj = getElement( parent, key );
      if ( !pobj )
        return;
      if ( !PyFloat_Check( pobj ) ) {
        Py_DECREF( pobj );
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      }
      out = PyFloat_AsDouble( pobj );
      Py_DECREF( pobj );
    }

    void
    PythonHandler::getParameter( PyObject* parent, const char* key, std::string& out )
    {
      PyObject* pobj = getElement( parent, key );
      if ( !pobj )
        return;
      if ( !PyString_Check( pobj ) ) {
        Py_DECREF( pobj );
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      }
      out = decode( pobj );
    }

    void
    PythonHandler::getParameter( PyObject* parent, const char* key, std::vector<std::string>& out )
    {
      PyObject* pobj = getElement( parent, key );
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) ) {
        Py_DECREF( pobj );
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      }
      for ( Py_ssize_t i = 0; i < PyTuple_Size( pobj ); ++i ) {
        PyObject* pit = PyTuple_GetItem( pobj, i );
        out.emplace_back( decode( pit ) );
      }
      Py_DECREF( pobj );
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
      return _PyInt_AsInt( obj );
#else
      return PyLong_AsLong( obj );
#endif
    }
  }
}

#endif

