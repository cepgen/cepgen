#include "ConfigReader.h"

#ifdef LIBCONFIG

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for CepGen input cards

    ConfigReader::ConfigReader( const char* file )
    {
      libconfig::Config cfg;
      try { cfg.readFile( file ); } catch ( const libconfig::ParseException& pe ) {
        FatalError( Form( "Failed to parse the configuration card \"%s\".\n\tParser error: %s (L:%d)", file, pe.getError(), pe.getLine() ) );
      }
      try {
        const libconfig::Setting& root = cfg.getRoot();
        const libconfig::Setting& proc = root["process"];

        //--- type of process to consider
        const std::string proc_name = proc["name"];
        if ( proc_name == "lpair" ) params_.setProcess( new Process::GamGamLL );
        else if ( proc_name == "pptoll" ) params_.setProcess( new Process::PPtoLL );
        else FatalError( Form( "Unrecognised process: %s", proc_name.c_str() ) );

        //--- process mode
        int int_mode; std::string str_mode;
        if ( proc.lookupValue( "mode", int_mode ) ) {
          params_.kinematics.mode = (Kinematics::ProcessMode)int_mode;
        }
        else if ( proc.lookupValue( "mode", str_mode ) ) {
          if ( str_mode == "elastic/elastic" ) params_.kinematics.mode = Kinematics::ProcessMode::ElasticElastic;
          else if ( str_mode == "elastic/inelastic" ) params_.kinematics.mode = Kinematics::ProcessMode::ElasticInelastic;
          else if ( str_mode == "inelastic/elastic" ) params_.kinematics.mode = Kinematics::ProcessMode::InelasticElastic;
          else if ( str_mode == "inelastic/inelastic" ) params_.kinematics.mode = Kinematics::ProcessMode::InelasticInelastic;
          else FatalError( Form( "Unrecognised interaction mode: %s", str_mode.c_str() ) );
        }

        //--- process kinematics
        if ( proc.exists( "in_kinematics" ) ) parseIncomingKinematics( proc["in_kinematics"] );
        if ( proc.exists( "out_kinematics" ) ) parseOutgoingKinematics( proc["out_kinematics"] );

        //--- generation parameters
        if ( root.exists( "vegas" ) ) parseVegas( root["vegas"] );
        if ( root.exists( "generator" ) ) parseGenerator( root["generator"] );

        //--- taming functions
        if ( proc.exists( "taming_functions" ) ) parseTamingFunctions( proc["taming_functions"] );

      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      } catch ( const libconfig::SettingTypeException& te ) {
        FatalError( Form( "Field \"%s\" has wrong type.", te.getPath() ) );
      }
    }

    void
    ConfigReader::parseIncomingKinematics( const libconfig::Setting& kin )
    {
      try {
        if ( kin.exists( "beam1_pz" ) ) params_.kinematics.in1p = (double)kin["beam1_pz"];
        if ( kin.exists( "beam2_pz" ) ) params_.kinematics.in2p = (double)kin["beam2_pz"];
        if ( kin.exists( "structure_functions" ) ) {
          std::string sf = kin["structure_functions" ];
          if ( sf == "electron" ) params_.remnant_mode = Electron;
          else if ( sf == "elastic proton" ) params_.remnant_mode = ElasticProton;
          else if ( sf == "Suri-Yennie" ) params_.remnant_mode = SuriYennie;
          else if ( sf == "Suri-Yennie;lowQ2" ) params_.remnant_mode = SuriYennieLowQ2;
          else if ( sf == "Szczurek-Uleshchenko" ) params_.remnant_mode = SzczurekUleshchenko;
          else if ( sf == "Fiore;valence" ) params_.remnant_mode = FioreVal;
          else if ( sf == "Fiore;sea" ) params_.remnant_mode = FioreSea;
          else if ( sf == "Fiore" ) params_.remnant_mode = Fiore;
          else FatalError( Form( "Invalid structure functions mode: %s", sf.c_str() ) );
        }
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      } catch ( const libconfig::SettingTypeException& te ) {
        FatalError( Form( "Field \"%s\" has wrong type.", te.getPath() ) );
      }
    }

    void
    ConfigReader::parseOutgoingKinematics( const libconfig::Setting& kin )
    {
      try {
        if ( kin.exists( "pair" ) ) params_.kinematics.pair = (Particle::ParticleCode)(int)kin["pair"];
        if ( kin.exists( "cuts_mode" ) ) params_.kinematics.cuts_mode = (Kinematics::Cuts)(int)kin["cuts_mode"];
        if ( kin.exists( "min_pt" ) ) params_.kinematics.pt_min = (double)kin["min_pt"];
        if ( kin.exists( "max_pt" ) ) params_.kinematics.pt_max = (double)kin["max_pt"];
        if ( kin.exists( "min_energy" ) ) params_.kinematics.e_min = (double)kin["min_energy"];
        if ( kin.exists( "max_energy" ) ) params_.kinematics.e_max = (double)kin["max_energy"];
        if ( kin.exists( "min_eta" ) ) params_.kinematics.eta_min = (double)kin["min_eta"];
        if ( kin.exists( "max_eta" ) ) params_.kinematics.eta_max = (double)kin["max_eta"];
        if ( kin.exists( "min_mx" ) ) params_.kinematics.mx_min = (double)kin["min_mx"];
        if ( kin.exists( "max_mx" ) ) params_.kinematics.mx_max = (double)kin["max_mx"];
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      } catch ( const libconfig::SettingTypeException& te ) {
        FatalError( Form( "Field \"%s\" has wrong type.", te.getPath() ) );
      }
    }

    void
    ConfigReader::parseVegas( const libconfig::Setting& veg )
    {
      try {
        if ( veg.exists( "num_points" ) ) params_.vegas.npoints = (int)veg["num_points"];
        if ( veg.exists( "num_integration_calls" ) ) params_.vegas.ncvg = (int)veg["num_integration_calls"];
        if ( veg.exists( "num_integration_iterations" ) ) params_.vegas.itvg = (int)veg["num_integration_iterations"];
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      }
    }

    void
    ConfigReader::parseGenerator( const libconfig::Setting& gen )
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
    ConfigReader::parseTamingFunctions( const libconfig::Setting& tf )
    {
      if ( !tf.isList() ) FatalError( "The taming functions definition must be wrapped within a list!" );
      for ( unsigned short i = 0; i < tf.getLength(); ++i ) {
        const std::string variable = tf[i]["variable"], expression = tf[i]["expression"];
        params_.taming_functions.emplace_back( variable, expression );
      }
    }

    void
    ConfigReader::writeProcess( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& proc = root.add( "process", libconfig::Setting::TypeGroup );
      proc.add( "name", libconfig::Setting::TypeString ) = params->processName();
      std::ostringstream os; os << params->kinematics.mode;
      proc.add( "mode", libconfig::Setting::TypeString ) = os.str();
    }

    void
    ConfigReader::writeIncomingKinematics( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& kin = root.add( "in_kinematics", libconfig::Setting::TypeGroup );
      kin.add( "beam1_pz", libconfig::Setting::TypeFloat ) = params->kinematics.in1p;
      kin.add( "beam2_pz", libconfig::Setting::TypeFloat ) = params->kinematics.in2p;
      std::ostringstream os; os << params->remnant_mode;
      kin.add( "structure_function", libconfig::Setting::TypeString ) = os.str();
    }

    void
    ConfigReader::writeOutgoingKinematics( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& kin = root.add( "out_kinematics", libconfig::Setting::TypeGroup );
      kin.add( "pair", libconfig::Setting::TypeInt ) = (int)params->kinematics.pair;
      kin.add( "cuts_mode", libconfig::Setting::TypeInt ) = (int)params->kinematics.cuts_mode;
      kin.add( "min_pt", libconfig::Setting::TypeFloat ) = params->kinematics.pt_min;
      kin.add( "max_pt", libconfig::Setting::TypeFloat ) = params->kinematics.pt_max;
      kin.add( "min_energy", libconfig::Setting::TypeFloat ) = params->kinematics.e_min;
      kin.add( "max_energy", libconfig::Setting::TypeFloat ) = params->kinematics.e_max;
      kin.add( "min_eta", libconfig::Setting::TypeFloat ) = params->kinematics.eta_min;
      kin.add( "max_eta", libconfig::Setting::TypeFloat ) = params->kinematics.eta_max;
      kin.add( "min_mx", libconfig::Setting::TypeFloat ) = params->kinematics.mx_min;
      kin.add( "max_mx", libconfig::Setting::TypeFloat ) = params->kinematics.mx_max;
    }

    void
    ConfigReader::writeTamingFunctions( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& tf = root.add( "taming_functions", libconfig::Setting::TypeList );
      for ( std::vector<Parameters::TamingFunction>::const_iterator it = params->taming_functions.begin(); it != params->taming_functions.end(); ++it ) {
        libconfig::Setting& fun = tf.add( libconfig::Setting::TypeGroup );
        fun.add( "variable", libconfig::Setting::TypeString ) = it->variable;
        fun.add( "expression", libconfig::Setting::TypeString ) = it->expression;
      }
    }

    void
    ConfigReader::writeVegas( const Parameters* params, libconfig::Setting& root )
    {
      libconfig::Setting& veg = root.add( "vegas", libconfig::Setting::TypeGroup );
      veg.add( "num_points", libconfig::Setting::TypeInt ) = (int)params->vegas.npoints;
      veg.add( "num_integration_calls", libconfig::Setting::TypeInt ) = (int)params->vegas.ncvg;
      veg.add( "num_integration_iterations", libconfig::Setting::TypeInt ) = (int)params->vegas.itvg;
    }

    void
    ConfigReader::writeGenerator( const Parameters* params, libconfig::Setting& root )
    {
      if ( !params->generation.enabled ) return;
      libconfig::Setting& gen = root.add( "generator", libconfig::Setting::TypeGroup );
      gen.add( "num_events", libconfig::Setting::TypeInt ) = (int)params->generation.maxgen;
      gen.add( "print_every", libconfig::Setting::TypeInt ) = (int)params->generation.gen_print_every;
    }

    void
    ConfigReader::store( const Parameters* params, const char* file )
    {
      libconfig::Config cfg;
      libconfig::Setting& root = cfg.getRoot();
      writeProcess( params, root );
      writeIncomingKinematics( params, root["process"] );
      writeOutgoingKinematics( params, root["process"] );
      writeTamingFunctions( params, root["process"] );
      writeVegas( params, root );
      writeGenerator( params, root );
      cfg.writeFile( file );
    }
  }
}

#endif
