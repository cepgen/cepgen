#include "ConfigHandler.h"
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
    ConfigHandler::ConfigHandler( const char* file )
    {
#ifdef LIBCONFIG
      libconfig::Config cfg;
      try { cfg.readFile( file ); } catch ( const libconfig::ParseException& pe ) {
        FatalError( Form( "Failed to parse the configuration card \"%s\".\n\tParser error: %s (at line %d)", file, pe.getError(), pe.getLine() ) );
      }
      try {
        const libconfig::Setting& root = cfg.getRoot();
        const libconfig::Setting& proc = root["process"];

        //--- type of process to consider
        const char* proc_name = proc["name"]; const std::string str_proc_name = proc_name;
        if ( str_proc_name == "lpair" ) params_.setProcess( new Process::GamGamLL );
        else if ( str_proc_name == "pptoll" ) params_.setProcess( new Process::PPtoLL );
        else if ( str_proc_name == "pptoww" ) params_.setProcess( new Process::PPtoWW );
        else FatalError( Form( "Unrecognised process: %s", proc_name ) );

        //--- process mode
        int int_mode; const char* mode;
        if ( proc.lookupValue( "mode", int_mode ) ) {
          params_.kinematics.mode = (Kinematics::ProcessMode)int_mode;
        }
        else if ( proc.lookupValue( "mode", mode ) ) {
          const std::string str_mode = mode;
          if ( str_mode == "elastic/elastic" ) params_.kinematics.mode = Kinematics::ProcessMode::ElasticElastic;
          else if ( str_mode == "elastic/inelastic" ) params_.kinematics.mode = Kinematics::ProcessMode::ElasticInelastic;
          else if ( str_mode == "inelastic/elastic" ) params_.kinematics.mode = Kinematics::ProcessMode::InelasticElastic;
          else if ( str_mode == "inelastic/inelastic" ) params_.kinematics.mode = Kinematics::ProcessMode::InelasticInelastic;
          else FatalError( Form( "Unrecognised interaction mode: %s", mode ) );
        }

        //--- process kinematics
        if ( proc.exists( "in_kinematics" ) ) parseIncomingKinematics( proc["in_kinematics"] );
        if ( proc.exists( "out_kinematics" ) ) parseOutgoingKinematics( proc["out_kinematics"] );

        //--- hadroniser parameters
        if ( root.exists( "hadroniser" ) ) parseHadroniser( root["hadroniser"] );

        //--- generation parameters
        if ( root.exists( "integrator" ) ) parseIntegrator( root["integrator"] );
        if ( root.exists( "vegas" ) ) parseIntegrator( root["vegas"] ); //FIXME for backward compatibility
        if ( root.exists( "generator" ) ) parseGenerator( root["generator"] );

        //--- taming functions
        if ( proc.exists( "taming_functions" ) ) parseTamingFunctions( proc["taming_functions"] );

      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      } catch ( const libconfig::SettingTypeException& te ) {
        FatalError( Form( "Field \"%s\" has wrong type.", te.getPath() ) );
      }
#else
      InWarning( "libconfig++ is not present on this machine" );
#endif
    }

#ifdef LIBCONFIG
    void
    ConfigHandler::parseIncomingKinematics( const libconfig::Setting& kin )
    {
      try {
        if ( kin.exists( "beam1_pz" ) ) params_.kinematics.inp.first = (double)kin["beam1_pz"];
        if ( kin.exists( "beam2_pz" ) ) params_.kinematics.inp.second = (double)kin["beam2_pz"];
        if ( kin.exists( "structure_functions" ) ) {
          const char* sf = kin["structure_functions" ]; const std::string sf_str = sf;
          if ( sf_str == "electron" ) params_.kinematics.structure_functions = StructureFunctions::Electron;
          else if ( sf_str == "elastic proton" ) params_.kinematics.structure_functions = StructureFunctions::ElasticProton;
          else if ( sf_str == "Suri-Yennie" ) params_.kinematics.structure_functions = StructureFunctions::SuriYennie;
          else if ( sf_str == "Szczurek-Uleshchenko" ) params_.kinematics.structure_functions = StructureFunctions::SzczurekUleshchenko;
          else if ( sf_str == "Fiore" ) params_.kinematics.structure_functions = StructureFunctions::FioreBrasse;
          else if ( sf_str == "ALLM" ) params_.kinematics.structure_functions = StructureFunctions::ALLM97;
          else if ( sf_str == "ALLM;91" ) params_.kinematics.structure_functions = StructureFunctions::ALLM91;
          else if ( sf_str == "ALLM;97" ) params_.kinematics.structure_functions = StructureFunctions::ALLM97;
          //else if ( sf_str == "ALLM;HHT" ) params_.kinematics.structure_functions = StructureFunctions::ALLM_HHT;
          //else if ( sf_str == "ALLM;HHT-FT" ) params_.kinematics.structure_functions = StructureFunctions::ALLM_HHT_FT;
          else if ( sf_str == "ALLM;GD07p" ) params_.kinematics.structure_functions = StructureFunctions::GD07p;
          else if ( sf_str == "ALLM;GD11p" ) params_.kinematics.structure_functions = StructureFunctions::GD11p;
          else if ( sf_str == "Schaefer" ) params_.kinematics.structure_functions = StructureFunctions::Schaefer;
          else FatalError( Form( "Invalid structure functions mode: %s", sf ) );
        }
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      } catch ( const libconfig::SettingTypeException& te ) {
        FatalError( Form( "Field \"%s\" has wrong type.", te.getPath() ) );
      }
    }

    void
    ConfigHandler::parseOutgoingKinematics( const libconfig::Setting& kin )
    {
      try {
        if ( kin.exists( "pair" ) ) {
          ParticleCode pair = (ParticleCode)(int)kin["pair"];
          params_.kinematics.central_system = { pair, pair };
        }
        if ( kin.exists( "min_pt" ) ) params_.kinematics.cuts.central[Cuts::pt_single].min() = (double)kin["min_pt"];
        if ( kin.exists( "max_pt" ) ) params_.kinematics.cuts.central[Cuts::pt_single].max() = (double)kin["max_pt"];
        if ( kin.exists( "min_ptdiff" ) ) params_.kinematics.cuts.central[Cuts::pt_diff].min() = (double)kin["min_ptdiff"];
        if ( kin.exists( "max_ptdiff" ) ) params_.kinematics.cuts.central[Cuts::pt_diff].max() = (double)kin["max_ptdiff"];
        if ( kin.exists( "min_rapiditydiff" ) ) params_.kinematics.cuts.central[Cuts::rapidity_diff].min() = (double)kin["min_rapiditydiff"];
        if ( kin.exists( "max_rapiditydiff" ) ) params_.kinematics.cuts.central[Cuts::rapidity_diff].max() = (double)kin["max_rapiditydiff"];
        if ( kin.exists( "min_energy" ) ) params_.kinematics.cuts.central[Cuts::energy_single].min() = (double)kin["min_energy"];
        if ( kin.exists( "max_energy" ) ) params_.kinematics.cuts.central[Cuts::energy_single].max() = (double)kin["max_energy"];
        if ( kin.exists( "min_eta" ) ) params_.kinematics.cuts.central[Cuts::eta_single].min() = (double)kin["min_eta"];
        if ( kin.exists( "max_eta" ) ) params_.kinematics.cuts.central[Cuts::eta_single].max() = (double)kin["max_eta"];
        if ( kin.exists( "min_rapidity" ) ) params_.kinematics.cuts.central[Cuts::rapidity_single].min() = (double)kin["min_rapidity"];
        if ( kin.exists( "max_rapidity" ) ) params_.kinematics.cuts.central[Cuts::rapidity_single].max() = (double)kin["max_rapidity"];
        if ( kin.exists( "min_mx" ) ) params_.kinematics.cuts.remnants[Cuts::mass].min() = (double)kin["min_mx"];
        if ( kin.exists( "max_mx" ) ) params_.kinematics.cuts.remnants[Cuts::mass].max() = (double)kin["max_mx"];
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      } catch ( const libconfig::SettingTypeException& te ) {
        FatalError( Form( "Field \"%s\" has wrong type.", te.getPath() ) );
      }
    }

    void
    ConfigHandler::parseIntegrator( const libconfig::Setting& integr )
    {
      try {
        if ( integr.exists( "algorithm" ) ) {
          const std::string algo = integr["algorithm"];
          if ( algo == "Vegas" ) params_.integrator.type = Integrator::Vegas;
          else if ( algo == "MISER" ) params_.integrator.type = Integrator::MISER;
        }
        if ( integr.exists( "num_points" ) ) params_.integrator.npoints = (int)integr["num_points"];
        if ( integr.exists( "num_integration_calls" ) ) params_.integrator.ncvg = (int)integr["num_integration_calls"];
        if ( integr.exists( "num_integration_iterations" ) ) params_.integrator.itvg = (int)integr["num_integration_iterations"];
        if ( integr.exists( "seed" ) ) params_.integrator.seed = (unsigned long)integr["seed"];
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      }
    }

    void
    ConfigHandler::parseGenerator( const libconfig::Setting& gen )
    {
      params_.generation.enabled = true;
      try {
        if ( gen.exists( "num_events" ) ) params_.generation.maxgen = (int)gen["num_events"];
        if ( gen.exists( "print_every" ) ) params_.generation.gen_print_every = (int)gen["print_every"];
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      }
    }

    void
    ConfigHandler::parseTamingFunctions( const libconfig::Setting& tf )
    {
      if ( !tf.isList() ) FatalError( "The taming functions definition must be wrapped within a list!" );
      for ( unsigned short i = 0; i < tf.getLength(); ++i ) {
        params_.taming_functions.add( tf[i]["variable"], tf[i]["expression"] );
      }
    }

    void
    ConfigHandler::parseHadroniser( const libconfig::Setting& hadr )
    {
      try {
        std::string hadr_name = hadr["name"];
        if ( hadr_name == "pythia8" ) {
#ifdef PYTHIA8
          Hadroniser::Pythia8Hadroniser* pythia8 = new Hadroniser::Pythia8Hadroniser;
          long long seed = ( hadr.exists( "seed" ) ) ? hadr["seed"] : -1ll;
          pythia8->setSeed( seed );
          pythia8->readString( Form( "Beams:idA = %d", params_.kinematics.inpdg.first ) );
          pythia8->readString( Form( "Beams:idB = %d", params_.kinematics.inpdg.second ) );
          pythia8->readString( Form( "Beams:eCM = %.2f", params_.kinematics.sqrtS() ) );
          if ( hadr.exists( "pythiaPreConfiguration" ) ) {
            libconfig::Setting& configs = hadr["pythiaPreConfiguration"];
            if ( !configs.isList() ) throw libconfig::SettingTypeException( configs );
            for ( unsigned short i = 0; i < configs.getLength(); ++i ) {
              std::string config = configs[i];
              pythia8->readString( config );
            }
          }
          pythia8->init();
          if ( hadr.exists( "pythiaConfiguration" ) ) {
            libconfig::Setting& configs = hadr["pythiaConfiguration"];
            if ( !configs.isList() ) throw libconfig::SettingTypeException( configs );
            for ( unsigned short i = 0; i < configs.getLength(); ++i ) {
              std::string config = configs[i];
              pythia8->readString( config );
            }
          }
          params_.setHadroniser( pythia8 );
#endif
        }
      } catch ( const libconfig::SettingTypeException& nfe ) {
        FatalError( Form( "Wrong type to the field \"%s\".", nfe.getPath() ) );
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      }
    }

    void
    ConfigHandler::writeProcess( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& proc = root.add( "process", libconfig::Setting::TypeGroup );
      proc.add( "name", libconfig::Setting::TypeString ) = params->processName();
      std::ostringstream os; os << params->kinematics.mode;
      proc.add( "mode", libconfig::Setting::TypeString ) = os.str();
    }

    void
    ConfigHandler::writeIncomingKinematics( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& kin = root.add( "in_kinematics", libconfig::Setting::TypeGroup );
      kin.add( "beam1_pz", libconfig::Setting::TypeFloat ) = params->kinematics.inp.first;
      kin.add( "beam2_pz", libconfig::Setting::TypeFloat ) = params->kinematics.inp.second;
      std::ostringstream os; os << params->kinematics.structure_functions;
      kin.add( "structure_function", libconfig::Setting::TypeString ) = os.str();
    }

    void
    ConfigHandler::writeOutgoingKinematics( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& kin = root.add( "out_kinematics", libconfig::Setting::TypeGroup );
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
    ConfigHandler::writeTamingFunctions( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& tf = root.add( "taming_functions", libconfig::Setting::TypeList );
      for ( std::map<std::string,TamingFunction>::const_iterator it = params->taming_functions.begin(); it != params->taming_functions.end(); ++it ) {
        libconfig::Setting& fun = tf.add( libconfig::Setting::TypeGroup );
        fun.add( "variable", libconfig::Setting::TypeString ) = it->first;
        fun.add( "expression", libconfig::Setting::TypeString ) = it->second.expression;
      }
    }

    void
    ConfigHandler::writeIntegrator( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& integr = root.add( "integrator", libconfig::Setting::TypeGroup );
      std::ostringstream os; os << params->integrator.type;
      integr.add( "algorithm", libconfig::Setting::TypeString ) = os.str();
      integr.add( "num_points", libconfig::Setting::TypeInt ) = (int)params->integrator.npoints;
      integr.add( "num_integration_calls", libconfig::Setting::TypeInt ) = (int)params->integrator.ncvg;
      integr.add( "num_integration_iterations", libconfig::Setting::TypeInt ) = (int)params->integrator.itvg;
      integr.add( "seed", libconfig::Setting::TypeInt64 ) = (long)params->integrator.seed;
    }

    void
    ConfigHandler::writeGenerator( const Parameters* params, libconfig::Setting& root )
    {
      if ( !params->generation.enabled ) return;
      libconfig::Setting& gen = root.add( "generator", libconfig::Setting::TypeGroup );
      gen.add( "num_events", libconfig::Setting::TypeInt ) = (int)params->generation.maxgen;
      gen.add( "print_every", libconfig::Setting::TypeInt ) = (int)params->generation.gen_print_every;
    }
#endif

    void
    ConfigHandler::store( const Parameters* params, const char* file )
    {
#ifdef LIBCONFIG
      libconfig::Config cfg;
      libconfig::Setting& root = cfg.getRoot();
      writeProcess( params, root );
      writeIncomingKinematics( params, root["process"] );
      writeOutgoingKinematics( params, root["process"] );
      writeTamingFunctions( params, root["process"] );
      writeIntegrator( params, root );
      writeGenerator( params, root );
      cfg.writeFile( file );
#endif
    }
  }
}

