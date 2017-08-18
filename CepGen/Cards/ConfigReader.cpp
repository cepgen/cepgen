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
          else if ( sf == "elastic-proton" ) params_.remnant_mode = ElasticProton;
          else if ( sf == "suri-yennie" ) params_.remnant_mode = SuriYennie;
          else if ( sf == "suri-yennie-lowQ2" ) params_.remnant_mode = SuriYennieLowQ2;
          else if ( sf == "szczurek-uleshchenko" ) params_.remnant_mode = SzczurekUleshchenko;
          else if ( sf == "fiore-valence" ) params_.remnant_mode = FioreVal;
          else if ( sf == "fiore-sea" ) params_.remnant_mode = FioreSea;
          else if ( sf == "fiore" ) params_.remnant_mode = Fiore;
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
    ConfigReader::store( const char* file ) const
    {
    }
  }
}

#endif
