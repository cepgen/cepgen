#include "CepGen/Cards/Handler.h"

namespace CepGen
{
  namespace Cards
  {
    //----- specialization for LPAIR input cards

    template<>
    LpairReader::Handler( const char* file )
    {
      std::ifstream f( file, std::fstream::in );
      if ( !f.is_open() ) {
        FatalError( Form( "Failed to parse file \"%s\"", file ) );
        return;
      }

      std::string key, value;
      const unsigned int wdth = 50;
      Debugging( Form( "File '%s' succesfully opened!", file ) );
      std::ostringstream os;
      os << "Configuration file content :" << "\n";

      os << std::left;
      while ( f >> key >> value ) {
        //os << std::setw( wdth ) << "[" << key << "] = " << value << std::endl;
        //if ( strncmp( key.c_str(), "#" ) == 0 ) continue; // FIXME need to ensure there is no extra space before !
        if ( key[0] == '#' ) continue;
        if ( key == "IEND" ) {
          int iend = (int)atoi( value.c_str() );
          if (iend>1) {
            params_.generation = true;
          }
        }
        else if ( key == "DEBG" ) {
          Logger::get().level = static_cast<Logger::LoggingLevel>( atoi( value.c_str() ) );
        }
        else if ( key == "NCVG" ) {
          params_.ncvg = (int)atoi( value.c_str() );
          os << std::setw( wdth ) << " * Number of function calls:" << params_.ncvg << "\n";
        }
        else if ( key == "NCSG" ) {
          params_.npoints = (int)atoi( value.c_str() );
          os << std::setw( wdth ) << " * Number of points to probe:" << params_.npoints << "\n";
        }
        else if ( key == "ITVG" ) {
          params_.itvg = (int)atoi( value.c_str() );
          os << std::setw( wdth ) << " * Number of Vegas iterations:" << params_.itvg << "\n";
        }
        else if ( key == "INPP" ) {
          params_.in1p = (double)atof( value.c_str() );
          os << std::setw( wdth ) << " * Momentum (1st primary particle):" << params_.in1p << " GeV/c\n";
        }
        else if ( key == "INPE" ) {
          params_.in2p = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Momentum (2nd primary particle):" << params_.in1p << " GeV/c\n";
        }
        else if ( key == "PROC" ) {
          if ( value == "lpair" )       params_.setProcess( new Process::GamGamLL() );
          else if ( value == "pptoll" ) params_.setProcess( new Process::PPtoLL() );
          std::ostringstream proc_name; proc_name << params_.process.get();
          os << std::setw( wdth ) << " * Process:" << boldify( proc_name.str() ) << "\n";
        }
        else if ( key == "HADR" ) {
#ifdef PYTHIA6
          if ( value == "pythia6" ) params_.setHadroniser( new Hadroniser::Pythia6Hadroniser );
#endif
#ifdef PYTHIA8
          if ( value == "pythia8" ) params_.setHadroniser( new Hadroniser::Pythia8Hadroniser );
#endif
          os << std::setw( wdth ) << " * Hadroniser:" << ( ( params_.hadroniser!=0 ) ? params_.hadroniser->name() : colourise( "*** no hadroniser ***", Colour::Red ) ) << "\n";
        }
        else if ( key == "MODE" ) {
          params_.process_mode = static_cast<Kinematics::ProcessMode>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * Subprocess' mode:" << static_cast<unsigned int>( params_.process_mode ) << " --> " << params_.process_mode << "\n";
        }
        else if ( key == "PMOD" or key == "EMOD" ) {
          params_.remnant_mode = static_cast<StructureFunctions>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * Outgoing primary particles' mode:" << static_cast<unsigned int>( params_.remnant_mode )
          	 << " --> " << params_.remnant_mode << "\n";
        }
        else if ( key == "PAIR" ) {
          params_.pair = static_cast<Particle::ParticleCode>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * Outgoing particles' PDG id:" << static_cast<unsigned int>( params_.pair )
             << " --> " << params_.pair << "\n";
        }
        else if ( key == "MCUT" ) {
          params_.mcut = static_cast<Kinematics::Cuts>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * Set of cuts to apply on final products:" << params_.mcut << "\n";
        }
        else if (key == "PTCT") {
          params_.minpt = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal transverse momentum (single central outgoing particle):" << params_.minpt << " GeV/c\n";
        }
        else if (key == "MSCT") {
          params_.minmass = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal central system mass:" << params_.minmass << " GeV/c**2\n";
        }
        else if (key == "ECUT") {
          params_.minenergy = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal energy (single central outgoing particle):" << params_.minenergy << " GeV\n";
        }
        else if (key == "NGEN") {
          params_.maxgen = static_cast<unsigned int>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * Number of events to generate:" << boldify( params_.maxgen ) << "\n";
        }
        else if (key == "THMN") {
          //params_.mintheta = atof( value.c_str() );
          //params_.setThetaRange( atof( value.c_str() ), 0. ); // FIXME FIXME
          os << std::setw( wdth ) << " * Minimal polar production angle for the central particles" << etaToTheta( params_.mineta ) << "\n";
        }
        else if (key == "THMX") {
          //params_.maxtheta = atof( value.c_str() );
          //params_.setThetaRange( 0., atof( value.c_str() ) ); //FIXME FIXME
          os << std::setw( wdth ) << " * Maximal polar production angle for the central particles" << etaToTheta( params_.maxeta ) << "\n";
        }
        else if (key == "ETMN") {
          params_.mineta = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal pseudo-rapidity (central outgoing particles):" << params_.mineta << "\n";
        }
        else if (key == "ETMX") {
          params_.maxeta = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Maximal pseudo-rapidity (central outgoing particles):" << params_.maxeta << "\n";
        }
        else if (key == "Q2MN") {
          params_.minq2 = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal Q^2 (exchanged parton):" << params_.minq2 << " GeV^2\n";
        }
        else if (key == "Q2MX") {
          params_.maxq2 = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Maximal Q^2 (exchanged parton):" << params_.maxq2 << " GeV^2\n";
        }
        else if (key == "MXMN") {
          params_.minmx = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal invariant mass of proton remnants:" << params_.minmx << " GeV/c^2\n";
        }
        else if (key == "MXMX") {
          params_.maxmx = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Maximal invariant mass of proton remnants:" << params_.maxmx << " GeV/c^2\n";
        }
        else if (key == "GPDF") {
          params_.gpdf = static_cast<unsigned int>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * GPDF:" << params_.gpdf << "\n";
        }
        else if (key == "SPDF") {
          params_.spdf = static_cast<unsigned int>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * SPDF:" << params_.spdf << "\n";
        }
        else if (key == "QPDF") {
          params_.qpdf = static_cast<unsigned int>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * QPDF:" << params_.qpdf << "\n";
        }
        else {
          InWarning( Form( "Unrecognized argument: [%s] = %s", key.c_str(), value.c_str() ) );
        }
      }
      f.close();

      Information( os.str() );
    }

    template<> void
    LpairReader::store( const char* file )
    {
      std::ofstream f( file, std::fstream::out | std::fstream::trunc );
      if ( !f.is_open() ) {
        InError( Form( "Failed to open file \"%s\" for writing", file ) );
        return;
      }
      // ...
      if ( params_.itvg > 0 ) f << "ITVG  " << params_.itvg << std::endl;
      if ( params_.minenergy !=- 1 ) f << "ECUT  " << params_.minenergy << std::endl;
      if ( params_.minenergy !=- 1 ) f << "PTCT  " << params_.minpt << std::endl;
      // ...
      f.close();
    }
  }
}

