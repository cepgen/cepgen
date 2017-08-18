#include "ConfigReader.h"

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for CepGen input cards

    ConfigReader::ConfigReader( const char* file )
    {
      libconfig::Config cfg;
      try { cfg.readFile( file ); } catch ( const libconfig::ParseException ) {
        FatalError( Form( "Failed to parse the configuration card \"%s\".", file ) );
      }
      try {
        const libconfig::Setting& root = cfg.getRoot();
        const libconfig::Setting& proc = root["process"];

        //--- type of process to consider
        const std::string proc_name = proc["name"];
        if ( proc_name == "lpair" ) params_.setProcess( new Process::GamGamLL );
        if ( proc_name == "pptoll" ) params_.setProcess( new Process::PPtoLL );

        //--- process kinematics
        parseKinematics( proc["kinematics"] );
      } catch ( const libconfig::SettingNotFoundException& nfe ) {
        FatalError( Form( "Failed to retrieve the field \"%s\".", nfe.getPath() ) );
      }
    }

    void
    ConfigReader::parseKinematics( const libconfig::Setting& kin )
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
