#include "PythonHandler.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"
#include "CepGen/Processes/PPtoWW.h"

#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for CepGen input cards
    PythonHandler::PythonHandler( const char* file )
    {
      setenv( "PYTHONPATH", ".", 1 );
      wchar_t* program = Py_DecodeLocale( "cepgen-python", nullptr );
      Py_SetProgramName( program );
      PyMem_RawFree( program );
      Py_Initialize();
      std::string filename = file;
      PyObject* fn = encode( filename.substr( 0, filename.find_last_of( "." ) ).c_str() );
      cfg_ = PyImport_Import( fn );
      Py_CLEAR( fn );
      if ( cfg_ == nullptr )
        throwPythonError( Form( "Failed to parse the configuration card %s", file ).c_str(), FatalError );

      PyObject* process = PyObject_GetAttrString( cfg_, "process" );

      //--- type of process to consider
      PyObject* pproc_name = PyDict_GetItem( process, encode( "name" ) );
      const char* proc_name = decode( pproc_name );
      Py_CLEAR( pproc_name );
      const std::string str_proc_name = proc_name;
      if ( str_proc_name == "lpair" ) params_.setProcess( new Process::GamGamLL );
      else if ( str_proc_name == "pptoll" ) params_.setProcess( new Process::PPtoLL );
      else if ( str_proc_name == "pptoww" ) params_.setProcess( new Process::PPtoWW );
      else FatalError( Form( "Unrecognised process: %s", proc_name ) );

      //--- process mode
      PyObject* pproc_mode = PyDict_GetItem( process, encode( "mode" ) );
      if ( pproc_mode ) {
        if ( PyLong_Check( pproc_mode ) ) {
          int int_mode = PyLong_AsLong( pproc_mode );
          params_.kinematics.mode = (Kinematics::ProcessMode)int_mode;
        }
        Py_CLEAR( pproc_mode );
      }

      //--- process kinematics
      PyObject* pin_kinematics = getElement( process, "inKinematics" );
      if ( pin_kinematics ) parseIncomingKinematics( pin_kinematics );
      Py_CLEAR( pin_kinematics );

      PyObject* pout_kinematics = getElement( process, "outKinematics" );
      if ( pout_kinematics ) parseOutgoingKinematics( pout_kinematics );
      Py_CLEAR( pout_kinematics );

      //--- hadroniser parameters
      PyObject* phad = PyObject_GetAttrString( cfg_, "hadroniser" );
      if ( phad ) parseHadroniser( phad );
      Py_CLEAR( phad );

      //--- generation parameters
      PyObject* pint = PyObject_GetAttrString( cfg_, "integrator" );
      if ( pint ) parseIntegrator( pint );
      Py_CLEAR( pint );

      PyObject* pgen = PyObject_GetAttrString( cfg_, "generator" );
      if ( pgen ) parseGenerator( pgen );
      Py_CLEAR( pgen );

      //--- taming functions
      PyObject* ptam = PyObject_GetAttrString( cfg_, "tamingFunctions" );
      if ( ptam ) parseTamingFunctions( ptam );
      Py_CLEAR( ptam );
    }

    PythonHandler::~PythonHandler()
    {
      Py_CLEAR( cfg_ );
      Py_Finalize();
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
        Py_CLEAR( ppz );
      }
      PyObject* psqrts = getElement( kin, "cmEnergy" );
      if ( psqrts && PyFloat_Check( psqrts ) ) {
        params_.kinematics.setSqrtS( PyFloat_AsDouble( psqrts ) );
        Py_CLEAR( psqrts );
      }
      PyObject* pstrfun = getElement( kin, "structureFunctions" );
      if ( pstrfun ) {
        if ( PyTuple_Check( pstrfun ) && PyTuple_Size( pstrfun ) > 0 ) {
          const std::string sf_str = decode( PyTuple_GetItem( pstrfun, 0 ) );
          std::string sf_var;
          if ( PyTuple_Size( pstrfun ) > 1 ) sf_var = decode( PyTuple_GetItem( pstrfun, 1 ) );
          if ( sf_str == "electron" )
            params_.kinematics.structure_functions = StructureFunctions::Electron;
          else if ( sf_str == "elastic proton" )
            params_.kinematics.structure_functions = StructureFunctions::ElasticProton;
          else if ( sf_str == "Suri-Yennie" )
            params_.kinematics.structure_functions = StructureFunctions::SuriYennie;
          else if ( sf_str == "Szczurek-Uleshchenko" )
            params_.kinematics.structure_functions = StructureFunctions::SzczurekUleshchenko;
          else if ( sf_str == "Fiore" )
            params_.kinematics.structure_functions = StructureFunctions::FioreBrasse;
          else if ( sf_str == "ALLM" ) {
            params_.kinematics.structure_functions = StructureFunctions::ALLM97;
            if ( sf_var == "91" )
              params_.kinematics.structure_functions = StructureFunctions::ALLM91;
            else if ( sf_var == "97" )
              params_.kinematics.structure_functions = StructureFunctions::ALLM97;
            /*else if ( sf_var == "HHT" )
              params_.kinematics.structure_functions = StructureFunctions::ALLM_HHT;
            else if ( sf_var == "HHT-FT" )
              params_.kinematics.structure_functions = StructureFunctions::ALLM_HHT_FT;*/
            else if ( sf_var == "GD07p" )
              params_.kinematics.structure_functions = StructureFunctions::GD07p;
            else if ( sf_var == "GD11p" )
              params_.kinematics.structure_functions = StructureFunctions::GD11p;
          }
          else if ( sf_str == "LUXlike" )
            params_.kinematics.structure_functions = StructureFunctions::Schaefer;
          else FatalError( Form( "Invalid structure functions mode: %s", sf_str.c_str() ) );
        }
        Py_CLEAR( pstrfun );
      }
    }

    void
    PythonHandler::parseOutgoingKinematics( PyObject* kin )
    {
      PyObject* ppair = getElement( kin, "pair" );
      if ( ppair && PyLong_Check( ppair ) ) {
        ParticleCode pair = (ParticleCode)PyLong_AsLong( ppair );
        params_.kinematics.central_system = { pair, pair };
        Py_CLEAR( ppair );
      }
      else if ( ppair && PyTuple_Check( ppair ) ) {
        if ( PyTuple_Size( ppair ) != 2 )
          FatalError( "Invalid value for in_kinematics.pair!" );
        ParticleCode pair1 = (ParticleCode)PyLong_AsLong( PyTuple_GetItem( ppair, 0 ) );
        ParticleCode pair2 = (ParticleCode)PyLong_AsLong( PyTuple_GetItem( ppair, 1 ) );
        params_.kinematics.central_system = { pair1, pair2 };
        Py_CLEAR( ppair );
      }

      PyObject* pcuts = getElement( kin, "cuts" );
      if ( pcuts && PyDict_Check( pcuts ) ) parseParticlesCuts( pcuts );

      getLimits( kin, "pt", params_.kinematics.cuts.central[Cuts::pt_single] );
      getLimits( kin, "rapidity", params_.kinematics.cuts.central[Cuts::rapidity_single] );
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
        ParticleCode pdg = (ParticleCode)PyLong_AsLong( pkey );
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
      PyObject* palgo = getElement( integr, "algorithm" );
      if ( !palgo || !PyUnicode_Check( palgo ) )
        throwPythonError( "Invalid integration algorithm!" );
      std::string algo = decode( palgo );
      if ( algo == "Vegas" )
        params_.integrator.type = Integrator::Vegas;
      else if ( algo == "MISER" )
        params_.integrator.type = Integrator::MISER;
      else
        throwPythonError( Form( "Invalid integration algorithm: %s", algo.c_str() ).c_str() );
      Py_CLEAR( palgo );

      PyObject* pnp = getElement( integr, "numPoints" );
      if ( pnp ) {
        if ( PyLong_AsLong( pnp ) )params_.integrator.npoints = PyLong_AsLong( pnp );
        Py_CLEAR( pnp );
      }
      PyObject* pnc = getElement( integr, "numIntegrationCalls" );
      if ( pnc ) {
        if ( PyLong_AsLong( pnc ) )params_.integrator.ncvg = PyLong_AsLong( pnc );
        Py_CLEAR( pnc );
      }
      PyObject* pnit = getElement( integr, "numIntegrationIterations" );
      if ( pnit ) {
        if ( PyLong_AsLong( pnit ) )params_.integrator.itvg = PyLong_AsLong( pnit );
        Py_CLEAR( pnit );
      }
      PyObject* psd = getElement( integr, "seed" );
      if ( psd ) {
        if ( PyLong_AsLong( psd ) )params_.integrator.seed = PyLong_AsUnsignedLong( psd );
        Py_CLEAR( psd );
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
        if ( PyLong_AsLong( pnev ) )params_.generation.maxgen = PyLong_AsLong( pnev );
        Py_CLEAR( pnev );
      }
      PyObject* ppev = getElement( gen, "printEvery" );
      if ( ppev ) {
        if ( PyLong_AsLong( ppev ) )params_.generation.gen_print_every = PyLong_AsLong( ppev );
        Py_CLEAR( ppev );
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
        Py_CLEAR( pvar );
        Py_CLEAR( pexpr );
      }
    }

    void
    PythonHandler::parseHadroniser( PyObject* hadr )
    {
      if ( !PyDict_Check( hadr ) )
        throwPythonError( "Hadroniser object should be a dictionary!" );

      PyObject* pname = getElement( hadr, "name" );
      if ( !pname )
        throwPythonError( "Hadroniser name is required!" );

      std::string hadr_name = decode( pname );
      Py_CLEAR( pname );

      if ( hadr_name == "pythia8" ) {
#ifdef PYTHIA8
        Hadroniser::Pythia8Hadroniser* pythia8 = new Hadroniser::Pythia8Hadroniser;
        PyObject* pseed = getElement( hadr, "seed" );
        long long seed = -1ll;
        if ( pseed ) {
          if ( PyLong_Check( pseed ) ) seed = PyLong_AsLongLong( pseed );
          Py_CLEAR( pseed );
        }
        pythia8->setSeed( seed );
        pythia8->readString( Form( "Beams:idA = %d", params_.kinematics.inpdg.first ) );
        pythia8->readString( Form( "Beams:idB = %d", params_.kinematics.inpdg.second ) );
        pythia8->readString( Form( "Beams:eCM = %.2f", params_.kinematics.sqrtS() ) );
        PyObject* ppc = getElement( hadr, "pythiaPreConfiguration" );
        if ( ppc ) {
          if ( PyTuple_Check( ppc ) ) {
            for ( Py_ssize_t i = 0; i < PyTuple_Size( ppc ); ++i ) {
              PyObject* pln = PyTuple_GetItem( ppc, i );
              std::string config = decode( pln );
              Py_CLEAR( pln );
              pythia8->readString( config );
            }
          }
          Py_CLEAR( ppc );
        }
        pythia8->init();
        ppc = getElement( hadr, "pythiaConfiguration" );
        if ( ppc ) {
          if ( PyTuple_Check( ppc ) ) {
            for ( Py_ssize_t i = 0; i < PyTuple_Size( ppc ); ++i ) {
              PyObject* pln = PyTuple_GetItem( ppc, i );
              std::string config = decode( pln );
              Py_CLEAR( pln );
              pythia8->readString( config );
            }
          }
          Py_CLEAR( ppc );
        }
        ppc = getElement( hadr, "pythiaProcessConfiguration" );
        if ( ppc ) {
          if ( PyTuple_Check( ppc ) ) {
            for ( Py_ssize_t i = 0; i < PyTuple_Size( ppc ); ++i ) {
              PyObject* pln = PyTuple_GetItem( ppc, i );
              std::string config = decode( pln );
              Py_CLEAR( pln );
              pythia8->readString( config );
            }
          }
          Py_CLEAR( ppc );
        }
        params_.setHadroniser( pythia8 );
#endif
      }
    }

    /*void
    PythonHandler::writeProcess( const Parameters* params, PyObject* root )
    {
      PyObject* proc = root.add( "process", libconfig::Setting::TypeGroup );
      proc.add( "name", libconfig::Setting::TypeString ) = params->processName();
      std::ostringstream os; os << params->kinematics.mode;
      proc.add( "mode", libconfig::Setting::TypeString ) = os.str();
    }

    void
    PythonHandler::writeIncomingKinematics( const Parameters* params, PyObject* root )
    {
      PyObject* kin = root.add( "in_kinematics", libconfig::Setting::TypeGroup );
      kin.add( "beam1_pz", libconfig::Setting::TypeFloat ) = params->kinematics.inp.first;
      kin.add( "beam2_pz", libconfig::Setting::TypeFloat ) = params->kinematics.inp.second;
      std::ostringstream os; os << params->kinematics.structure_functions;
      kin.add( "structure_function", libconfig::Setting::TypeString ) = os.str();
    }

    void
    PythonHandler::writeOutgoingKinematics( const Parameters* params, PyObject* root )
    {
      PyObject* kin = root.add( "out_kinematics", libconfig::Setting::TypeGroup );
      if ( params->kinematics.central_system.size() > 0 )
        kin.add( "pair", libconfig::Setting::TypeInt ) = (int)params->kinematics.central_system[0];
      if ( params->kinematics.cuts.central.count( Cuts::pt_single ) ) {
        kin.add( "min_pt", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::pt_single ).min();
        kin.add( "max_pt", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::pt_single ).max();
      }
      if ( params->kinematics.cuts.central.count( Cuts::pt_diff ) ) {
        kin.add( "min_ptdiff", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::pt_diff ).min();
        kin.add( "max_ptdiff", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::pt_diff ).max();
      }
      if ( params->kinematics.cuts.central.count( Cuts::rapidity_diff ) ) {
        kin.add( "min_rapiditydiff", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::rapidity_diff ).min();
        kin.add( "max_rapiditydiff", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::rapidity_diff ).max();
      }
      if ( params->kinematics.cuts.central.count( Cuts::energy_single ) ) {
        kin.add( "min_energy", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::energy_single ).min();
        kin.add( "max_energy", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::energy_single ).max();
      }
      if ( params->kinematics.cuts.central.count( Cuts::eta_single ) ) {
        kin.add( "min_eta", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::eta_single ).min();
        kin.add( "max_eta", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.central.at( Cuts::eta_single ).max();
      }
      if ( params->kinematics.cuts.remnants.count( Cuts::mass ) ) {
        kin.add( "min_mx", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.remnants.at( Cuts::mass ).min();
        kin.add( "max_mx", libconfig::Setting::TypeFloat ) = params->kinematics.cuts.remnants.at( Cuts::mass ).max();
      }
    }

    void
    PythonHandler::writeTamingFunctions( const Parameters* params, PyObject* root )
    {
      PyObject* tf = root.add( "taming_functions", libconfig::Setting::TypeList );
      for ( std::map<std::string,TamingFunction>::const_iterator it = params->taming_functions.begin(); it != params->taming_functions.end(); ++it ) {
        PyObject* fun = tf.add( libconfig::Setting::TypeGroup );
        fun.add( "variable", libconfig::Setting::TypeString ) = it->first;
        fun.add( "expression", libconfig::Setting::TypeString ) = it->second.expression;
      }
    }

    void
    PythonHandler::writeIntegrator( const Parameters* params, PyObject* root )
    {
      PyObject* integr = root.add( "integrator", libconfig::Setting::TypeGroup );
      std::ostringstream os; os << params->integrator.type;
      integr.add( "algorithm", libconfig::Setting::TypeString ) = os.str();
      integr.add( "num_points", libconfig::Setting::TypeInt ) = (int)params->integrator.npoints;
      integr.add( "num_integration_calls", libconfig::Setting::TypeInt ) = (int)params->integrator.ncvg;
      integr.add( "num_integration_iterations", libconfig::Setting::TypeInt ) = (int)params->integrator.itvg;
      integr.add( "seed", libconfig::Setting::TypeInt64 ) = (long)params->integrator.seed;
    }

    void
    PythonHandler::writeGenerator( const Parameters* params, PyObject* root )
    {
      if ( !params->generation.enabled ) return;
      PyObject* gen = root.add( "generator", libconfig::Setting::TypeGroup );
      gen.add( "num_events", libconfig::Setting::TypeInt ) = (int)params->generation.maxgen;
      gen.add( "print_every", libconfig::Setting::TypeInt ) = (int)params->generation.gen_print_every;
    }

    void
    PythonHandler::store( const Parameters* params, const char* file )
    {
      libconfig::Config cfg;
      PyObject* root = cfg.getRoot();
      writeProcess( params, root );
      writeIncomingKinematics( params, root["process"] );
      writeOutgoingKinematics( params, root["process"] );
      writeTamingFunctions( params, root["process"] );
      writeIntegrator( params, root );
      writeGenerator( params, root );
      cfg.writeFile( file );
    }*/
    void
    PythonHandler::throwPythonError( const char* message, const ExceptionType& type ) const
    {
      PyObject* ptype = nullptr, *pvalue = nullptr, *ptraceback = nullptr;
      PyErr_Fetch( &ptype, &pvalue, &ptraceback );
      PyErr_Clear();
      PyErr_NormalizeException( &ptype, &pvalue, &ptraceback );
      std::ostringstream oss; oss << message;
      if ( ptype == nullptr ) {
        Py_DECREF( cfg_ );
        Py_Finalize();
        throw Exception( __PRETTY_FUNCTION__, oss.str().c_str(), type );
      }

      oss << "\n\tError: " << decode( PyObject_Str( pvalue ) );
      Py_DECREF( cfg_ );
      Py_Finalize();
      throw Exception( __PRETTY_FUNCTION__, oss.str().c_str(), type );
    }

    const char*
    PythonHandler::decode( PyObject* obj )
    {
      return _PyUnicode_AsString( obj );
    }

    PyObject*
    PythonHandler::encode( const char* str )
    {
      return PyUnicode_FromString( str );
    }

    PyObject*
    PythonHandler::getElement( PyObject* obj, const char* key )
    {
      PyObject* pout = nullptr;
      PyObject* nink = encode( key );
      pout = PyDict_GetItem( obj, nink );
      Py_CLEAR( nink );
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
        Py_CLEAR( pobj );
      }
    }
  }
}
