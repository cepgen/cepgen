#include "PythonHandler.h"
#include "CepGen/Core/Exception.h"

#ifdef PYTHON

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"
#include "CepGen/Processes/PPtoWW.h"

#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

#include <algorithm>

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
    std::string config = CepGen::Cards::PythonHandler::decode( pln );
    Py_DECREF( pln );
    py8->readString( config );
  }
  Py_DECREF( ppc );
}
#endif

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for CepGen input cards
    PythonHandler::PythonHandler( const char* file )
    {
      setenv( "PYTHONPATH", ".:..", 1 );
      std::string filename = getPythonPath( file );
      const size_t fn_len = filename.length()+1;
#ifdef PYTHON2
      {
        char* n_filename = new char[fn_len];
        sprintf( n_filename, "%s", filename.c_str() );
        Py_SetProgramName( n_filename );
        delete [] n_filename;
      }
#else
      {
        wchar_t* w_filename = new wchar_t[fn_len];
        swprintf( w_filename, fn_len, L"%s", filename.c_str() );
        Py_SetProgramName( w_filename );
        delete [] w_filename;
      }
#endif

      Py_InitializeEx( 0 );
      if ( !Py_IsInitialized() )
        throw Exception( __PRETTY_FUNCTION__, "Failed to initialise the python parser!", FatalError );

      Debugging( Form( "Initialised the Python cards parser\n\t"
                       "Python version: %s\n\t"
                       "Platform: %s", Py_GetVersion(), Py_GetPlatform() ) );

      PyObject* fn = encode( filename.c_str() );
      if ( !fn )
        throwPythonError( Form( "Failed to encode the configuration filename %s", filename.c_str() ).c_str() );

      PyObject* cfg = PyImport_Import( fn );
      Py_DECREF( fn );
      if ( !cfg )
        throwPythonError( Form( "Failed to parse the configuration card %s", file ).c_str(), FatalError );

      PyObject* process = PyObject_GetAttrString( cfg, "process" );

      //--- type of process to consider
      PyObject* pproc_name = PyDict_GetItem( process, encode( module_name_ ) );
      const char* proc_name = decode( pproc_name );
      Py_DECREF( pproc_name );
      const std::string str_proc_name = proc_name;
      if      ( str_proc_name == "lpair" )  params_.setProcess( new Process::GamGamLL );
      else if ( str_proc_name == "pptoll" ) params_.setProcess( new Process::PPtoLL );
      else if ( str_proc_name == "pptoww" ) params_.setProcess( new Process::PPtoWW );
      else FatalError( Form( "Unrecognised process: %s", proc_name ) );

      //--- process mode
      PyObject* pproc_mode = getElement( process, "mode" );
      if ( pproc_mode ) {
        if ( isInteger( pproc_mode ) ) {
          int int_mode = asInteger( pproc_mode );
          params_.kinematics.mode = (Kinematics::ProcessMode)int_mode;
        }
        Py_DECREF( pproc_mode );
      }

      //--- process kinematics
      PyObject* pin_kinematics = getElement( process, "inKinematics" );
      if ( pin_kinematics ) parseIncomingKinematics( pin_kinematics );
      Py_DECREF( pin_kinematics );

      PyObject* pout_kinematics = getElement( process, "outKinematics" );
      if ( pout_kinematics ) parseOutgoingKinematics( pout_kinematics );
      Py_DECREF( pout_kinematics );

      Py_DECREF( process );

      //--- hadroniser parameters
      PyObject* phad = PyObject_GetAttrString( cfg, "hadroniser" );
      if ( phad ) {
        parseHadroniser( phad );
        Py_DECREF( phad );
      }

      //--- generation parameters
      PyObject* pint = PyObject_GetAttrString( cfg, "integrator" );
      if ( pint ) {
        parseIntegrator( pint );
        Py_DECREF( pint );
      }

      PyObject* pgen = PyObject_GetAttrString( cfg, "generator" );
      if ( pgen ) {
        parseGenerator( pgen );
        Py_DECREF( pgen );
      }

      //--- taming functions
      PyObject* ptam = PyObject_GetAttrString( cfg, "tamingFunctions" );
      if ( ptam ) {
        parseTamingFunctions( ptam );
        Py_DECREF( ptam );
      }

      Py_DECREF( cfg );
#ifdef PYTHON2
      Py_Finalize();
#else
      if ( Py_IsInitialized() )
        throw Exception( __PRETTY_FUNCTION__, "Failed to unregister the python parser!", FatalError );
      if ( Py_FinalizeEx() != 0 )
        throw Exception( __PRETTY_FUNCTION__, "Failed to unregister the python parser!", FatalError );
#endif
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
      PyObject* psqrts = getElement( kin, "cmEnergy" );
      if ( psqrts ) {
        if ( PyFloat_Check( psqrts ) )
          params_.kinematics.setSqrtS( PyFloat_AsDouble( psqrts ) );
        Py_DECREF( psqrts );
      }
      PyObject* pstrfun = getElement( kin, "structureFunctions" );
      if ( pstrfun ) {
        if ( isInteger( pstrfun ) )
          params_.kinematics.structure_functions = (StructureFunctions::Type)asInteger( pstrfun );
        Py_DECREF( pstrfun );
      }
    }

    void
    PythonHandler::parseOutgoingKinematics( PyObject* kin )
    {
      PyObject* ppair = getElement( kin, "pair" );
      if ( ppair ) {
        if ( isInteger( ppair ) ) {
          ParticleCode pair = (ParticleCode)asInteger( ppair );
          params_.kinematics.central_system = { pair, pair };
        }
        Py_DECREF( ppair );
      }
      else if ( ppair && PyTuple_Check( ppair ) ) {
        if ( PyTuple_Size( ppair ) != 2 )
          FatalError( "Invalid value for in_kinematics.pair!" );
        ParticleCode pair1 = (ParticleCode)asInteger( PyTuple_GetItem( ppair, 0 ) );
        ParticleCode pair2 = (ParticleCode)asInteger( PyTuple_GetItem( ppair, 1 ) );
        params_.kinematics.central_system = { pair1, pair2 };
        Py_DECREF( ppair );
      }

      PyObject* pcuts = getElement( kin, "cuts" );
      if ( pcuts && PyDict_Check( pcuts ) ) parseParticlesCuts( pcuts );

      // for LPAIR/collinear matrix elements
      getLimits( kin, "q2", params_.kinematics.cuts.initial[Cuts::q2] );

      // for the kT factorised matrix elements
      getLimits( kin, "qt", params_.kinematics.cuts.initial[Cuts::qt] );
      getLimits( kin, "ptdiff", params_.kinematics.cuts.central[Cuts::pt_diff] );
      getLimits( kin, "rapiditydiff", params_.kinematics.cuts.central[Cuts::rapidity_diff] );

      getLimits( kin, "mx", params_.kinematics.cuts.remnants[Cuts::mass] );
    }

    void
    PythonHandler::parseParticlesCuts( PyObject* cuts )
    {
      PyObject* pkey = nullptr, *pvalue = nullptr;
      Py_ssize_t pos = 0;
      while ( PyDict_Next( cuts, &pos, &pkey, &pvalue ) ) {
        ParticleCode pdg = (ParticleCode)asInteger( pkey );
        getLimits( pvalue, "pt", params_.kinematics.cuts.central_particles[pdg][Cuts::pt_single] );
        getLimits( pvalue, "energy", params_.kinematics.cuts.central_particles[pdg][Cuts::energy_single] );
        getLimits( pvalue, "eta", params_.kinematics.cuts.central_particles[pdg][Cuts::eta_single] );
        getLimits( pvalue, "rapidity", params_.kinematics.cuts.central_particles[pdg][Cuts::rapidity_single] );
      }
    }

    void
    PythonHandler::parseIntegrator( PyObject* integr )
    {
      if ( !PyDict_Check( integr ) )
        throwPythonError( "Integrator object should be a dictionary!" );
      PyObject* palgo = getElement( integr, module_name_ );
      if ( !palgo )
        throwPythonError( "Invalid integration algorithm!" );
      std::string algo = decode( palgo );
      if      ( algo == "Plain" ) params_.integrator.type = Integrator::Plain;
      else if ( algo == "Vegas" ) {
        params_.integrator.type = Integrator::Vegas;
        PyObject* palpha = getElement( integr, "alpha" );
        if ( palpha ) {
          if ( PyFloat_Check( palpha ) )
            params_.integrator.vegas.alpha = PyFloat_AsDouble( palpha );
          Py_DECREF( palpha );
        }
        PyObject* piter = getElement( integr, "iterations" );
        if ( piter ) {
          if ( isInteger( piter ) ) {
            params_.integrator.vegas.iterations = asInteger( piter );
        }
          Py_DECREF( piter );
        }
        PyObject* pmode = getElement( integr, "mode" );
        if ( pmode ) {
          if ( isInteger( pmode ) )
            params_.integrator.vegas.mode = asInteger( pmode );
          Py_DECREF( pmode );
        }
        PyObject* pverb = getElement( integr, "verbosity" );
        if ( pverb ) {
          if ( isInteger( pverb ) )
            params_.integrator.vegas.verbose = asInteger( pverb );
          Py_DECREF( pverb );
        }
      }
      else if ( algo == "MISER" ) {
        params_.integrator.type = Integrator::MISER;
        PyObject* pestfrac = getElement( integr, "estimateFraction" );
        if ( pestfrac ) {
          if ( PyFloat_Check( pestfrac ) )
            params_.integrator.miser.estimate_frac = PyFloat_AsDouble( pestfrac );
          Py_DECREF( pestfrac );
        }
        PyObject* pmc = getElement( integr, "minCalls" );
        if ( pmc ) {
          if ( isInteger( pmc ) )
            params_.integrator.miser.min_calls = asInteger( pmc );
          Py_DECREF( pmc );
        }
        PyObject* pmcpb = getElement( integr, "minCallsPerBisection" );
        if ( pmcpb ) {
          if ( isInteger( pmcpb ) )
            params_.integrator.miser.min_calls_per_bisection = asInteger( pmcpb );
          Py_DECREF( pmcpb );
        }
        PyObject* palpha = getElement( integr, "alpha" );
        if ( palpha ) {
          if ( PyFloat_Check( palpha ) )
            params_.integrator.miser.alpha = PyFloat_AsDouble( palpha );
          Py_DECREF( palpha );
        }
        PyObject* pdither = getElement( integr, "dither" );
        if ( pdither ) {
          if ( PyFloat_Check( pdither ) )
            params_.integrator.miser.dither = PyFloat_AsDouble( pdither );
          Py_DECREF( pdither );
        }
      }
      else
        throwPythonError( Form( "Invalid integration algorithm: %s", algo.c_str() ).c_str() );
      Py_DECREF( palgo );

      PyObject* pnp = getElement( integr, "numPoints" );
      if ( pnp ) {
        if ( isInteger( pnp ) )
          params_.integrator.npoints = asInteger( pnp );
        Py_DECREF( pnp );
      }
      PyObject* pnc = getElement( integr, "numFunctionCalls" );
      if ( pnc ) {
        if ( isInteger( pnc ) )
          params_.integrator.ncvg = asInteger( pnc );
        Py_DECREF( pnc );
      }
      PyObject* psd = getElement( integr, "seed" );
      if ( psd ) {
        if ( isInteger( psd ) )
          params_.integrator.seed = PyLong_AsUnsignedLong( psd );
        Py_DECREF( psd );
      }
    }

    void
    PythonHandler::parseGenerator( PyObject* gen )
    {
      if ( !PyDict_Check( gen ) )
        throwPythonError( "Generation information object should be a dictionary!" );
      params_.generation.enabled = true;
      PyObject* pnev = getElement( gen, "numEvents" );
      if ( pnev ) {
        if ( isInteger( pnev ) )
          params_.generation.maxgen = asInteger( pnev );
        Py_DECREF( pnev );
      }
      PyObject* ppev = getElement( gen, "printEvery" );
      if ( ppev ) {
        if ( isInteger( ppev ) )
          params_.generation.gen_print_every = asInteger( ppev );
        Py_DECREF( ppev );
      }
    }

    void
    PythonHandler::parseTamingFunctions( PyObject* tf )
    {
      if ( !PyList_Check( tf ) )
        throwPythonError( "Taming functions list should be a list!" );

      for ( Py_ssize_t i = 0; i < PyList_Size( tf ); ++i ) {
        PyObject* pit = PyList_GetItem( tf, i );
        if ( !PyDict_Check( pit ) )
          throwPythonError( Form( "Item %d is invalid", i ).c_str() );
        PyObject* pvar = getElement( pit, "variable" ), *pexpr = getElement( pit, "expression" );
        params_.taming_functions.add( decode( pvar ), decode( pexpr ) );
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
      Py_DECREF( pname );

      if ( hadr_name == "pythia8" ) {
#ifdef PYTHIA8
        Hadroniser::Pythia8Hadroniser* pythia8 = new Hadroniser::Pythia8Hadroniser( params_ );
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
        params_.setHadroniser( pythia8 );
#else
        InWarning( "Pythia8 is not linked to this instance... "
                   "Ignoring this part of the configuration file." )
#endif
      }
    }

    //------------------------------------------------------------------

    std::string
    PythonHandler::getPythonPath( const char* file )
    {
      std::string s_filename = file;
      s_filename = s_filename.substr( 0, s_filename.find_last_of( "." ) );
      std::replace( s_filename.begin(), s_filename.end(), '/', '.' );
      return s_filename;
    }

    void
    PythonHandler::throwPythonError( const char* message, const ExceptionType& type )
    {
      PyObject* ptype = nullptr, *pvalue = nullptr, *ptraceback = nullptr;
      PyErr_Fetch( &ptype, &pvalue, &ptraceback );
      PyErr_Clear();
      PyErr_NormalizeException( &ptype, &pvalue, &ptraceback );
      std::ostringstream oss; oss << message;
      if ( ptype == nullptr ) {
        Py_Finalize();
        throw Exception( __PRETTY_FUNCTION__, oss.str().c_str(), type );
      }

      oss << "\n\tError: "
#ifdef PYTHON2
          << PyString_AsString( PyObject_Str( pvalue ) ); // deprecated in python v3+
#else
          << _PyUnicode_AsString( PyObject_Str( pvalue ) );
#endif
      Py_Finalize();
      throw Exception( __PRETTY_FUNCTION__, oss.str().c_str(), type );
    }

    const char*
    PythonHandler::decode( PyObject* obj )
    {
#ifdef PYTHON2
      const char* str = PyString_AsString( obj ); // deprecated in python v3+
#else
      const char* str = _PyUnicode_AsString( obj );
#endif
      if ( !str )
        throwPythonError( "Failed to decode a Python object!" );
      return str;
    }

    PyObject*
    PythonHandler::encode( const char* str )
    {
      PyObject* obj = PyUnicode_FromString( str );
      if ( !obj )
        throwPythonError( Form( "Failed to encode the following string:\n\t%s", str ).c_str() );
      return obj;
    }

    PyObject*
    PythonHandler::getElement( PyObject* obj, const char* key )
    {
      PyObject* nink = encode( key );
      PyObject* pout = PyDict_GetItem( obj, nink );
      Py_DECREF( nink );
      return pout;
    }

    void
    PythonHandler::getLimits( PyObject* obj, const char* key, Kinematics::Limits& lim )
    {
      PyObject* pobj = getElement( obj, key );
      if ( pobj ) {
        if ( !PyTuple_Check( pobj ) )
          FatalError( Form( "Invalid value retrieved for %s", key ) );
        if ( PyTuple_Size( pobj ) < 1 )
          FatalError( Form( "Invalid number of values unpacked for %s!", key ) );
        double min = PyFloat_AsDouble( PyTuple_GetItem( pobj, 0 ) );
        lim.min() = min;
        if ( PyTuple_Size( pobj ) > 1 ) {
          double max = PyFloat_AsDouble( PyTuple_GetItem( pobj, 1 ) );
          if ( max != -1 )
            lim.max() = max;
        }
        Py_DECREF( pobj );
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
      return _PyInt_AsInt( obj );
#else
      return PyLong_AsLong( obj );
#endif
    }
  }
}

#endif
