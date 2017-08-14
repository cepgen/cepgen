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
          params_.vegas.ncvg = (int)atoi( value.c_str() );
          os << std::setw( wdth ) << " * Number of function calls:" << params_.vegas.ncvg << "\n";
        }
        else if ( key == "NCSG" ) {
          params_.vegas.npoints = (int)atoi( value.c_str() );
          os << std::setw( wdth ) << " * Number of points to probe:" << params_.vegas.npoints << "\n";
        }
        else if ( key == "ITVG" ) {
          params_.vegas.itvg = (int)atoi( value.c_str() );
          os << std::setw( wdth ) << " * Number of Vegas iterations:" << params_.vegas.itvg << "\n";
        }
        else if ( key == "INPP" ) {
          params_.kinematics.in1p = (double)atof( value.c_str() );
          os << std::setw( wdth ) << " * Momentum (1st primary particle):" << params_.kinematics.in1p << " GeV/c\n";
        }
        else if ( key == "INPE" ) {
          params_.kinematics.in2p = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Momentum (2nd primary particle):" << params_.kinematics.in1p << " GeV/c\n";
        }
        else if ( key == "PROC" ) {
          if ( value == "lpair" )       params_.setProcess( new Process::GamGamLL() );
          else if ( value == "pptoll" ) params_.setProcess( new Process::PPtoLL() );
          std::ostringstream proc_name; proc_name << params_.process();
          os << std::setw( wdth ) << " * Process:" << boldify( proc_name.str() ) << "\n";
        }
        else if ( key == "HADR" ) {
#ifdef PYTHIA6
          if ( value == "pythia6" ) params_.setHadroniser( new Hadroniser::Pythia6Hadroniser );
#endif
#ifdef PYTHIA8
          if ( value == "pythia8" ) params_.setHadroniser( new Hadroniser::Pythia8Hadroniser );
#endif
          os << std::setw( wdth ) << " * Hadroniser:" << ( ( params_.hadroniser() != 0 ) ? params_.hadroniser()->name() : colourise( "*** no hadroniser ***", Colour::Red ) ) << "\n";
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
          params_.kinematics.pair = static_cast<Particle::ParticleCode>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * Outgoing particles' PDG id:" << static_cast<unsigned int>( params_.kinematics.pair )
             << " --> " << params_.kinematics.pair << "\n";
        }
        else if ( key == "MCUT" ) {
          params_.kinematics.cuts_mode = static_cast<Kinematics::Cuts>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * Set of cuts to apply on final products:" << params_.kinematics.cuts_mode << "\n";
        }
        else if ( key == "PTCT" ) {
          params_.kinematics.pt_min = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal transverse momentum (single central outgoing particle):" << params_.kinematics.pt_min << " GeV/c\n";
        }
        else if ( key == "MSCT" ) {
          params_.kinematics.mass_min = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal central system mass:" << params_.kinematics.mass_min << " GeV/c**2\n";
        }
        else if ( key == "ECUT" ) {
          params_.kinematics.e_min = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal energy (single central outgoing particle):" << params_.kinematics.e_min << " GeV\n";
        }
        else if ( key == "NGEN" ) {
          params_.maxgen = static_cast<unsigned int>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * Number of events to generate:" << boldify( params_.maxgen ) << "\n";
        }
        else if ( key == "THMN" ) {
          //params_.mintheta = atof( value.c_str() );
          //params_.setThetaRange( atof( value.c_str() ), 0. ); // FIXME FIXME
          os << std::setw( wdth ) << " * Minimal polar production angle for the central particles" << etaToTheta( params_.kinematics.eta_min ) << "\n";
        }
        else if ( key == "THMX" ) {
          //params_.maxtheta = atof( value.c_str() );
          //params_.setThetaRange( 0., atof( value.c_str() ) ); //FIXME FIXME
          os << std::setw( wdth ) << " * Maximal polar production angle for the central particles" << etaToTheta( params_.kinematics.eta_max ) << "\n";
        }
        else if ( key == "ETMN" ) {
          params_.kinematics.eta_min = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal pseudo-rapidity (central outgoing particles):" << params_.kinematics.eta_min << "\n";
        }
        else if ( key == "ETMX" ) {
          params_.kinematics.eta_max = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Maximal pseudo-rapidity (central outgoing particles):" << params_.kinematics.eta_max << "\n";
        }
        else if ( key == "Q2MN" ) {
          params_.kinematics.q2_min = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal Q^2 (exchanged parton):" << params_.kinematics.q2_min << " GeV^2\n";
        }
        else if ( key == "Q2MX" ) {
          params_.kinematics.q2_max = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Maximal Q^2 (exchanged parton):" << params_.kinematics.q2_max << " GeV^2\n";
        }
        else if ( key == "MXMN" ) {
          params_.kinematics.mx_min = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Minimal invariant mass of proton remnants:" << params_.kinematics.mx_min << " GeV/c^2\n";
        }
        else if ( key == "MXMX" ) {
          params_.kinematics.mx_max = static_cast<float>( atof( value.c_str() ) );
          os << std::setw( wdth ) << " * Maximal invariant mass of proton remnants:" << params_.kinematics.mx_max << " GeV/c^2\n";
        }
        else if ( key == "GPDF" ) {
          params_.pdflib.gpdf = static_cast<unsigned int>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * GPDF:" << params_.pdflib.gpdf << "\n";
        }
        else if ( key == "SPDF" ) {
          params_.pdflib.spdf = static_cast<unsigned int>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * SPDF:" << params_.pdflib.spdf << "\n";
        }
        else if ( key == "QPDF" ) {
          params_.pdflib.qpdf = static_cast<unsigned int>( atoi( value.c_str() ) );
          os << std::setw( wdth ) << " * QPDF:" << params_.pdflib.qpdf << "\n";
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
      if ( params_.vegas.itvg > 0 ) f << "ITVG  " << params_.vegas.itvg << std::endl;
      if ( params_.kinematics.e_min !=- 1 ) f << "ECUT  " << params_.kinematics.e_min << std::endl;
      if ( params_.kinematics.pt_min !=- 1 ) f << "PTCT  " << params_.kinematics.pt_min << std::endl;
      // ...
      f.close();
    }
  }
}

