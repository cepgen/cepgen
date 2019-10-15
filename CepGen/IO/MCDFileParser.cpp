#include "CepGen/IO/MCDFileParser.h"
#include "CepGen/Core/Exception.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace pdg
{
  void
  MCDFileParser::parse( const char* path )
  {
    std::ifstream ifile( path );
    std::string line;
    while ( std::getline( ifile, line ) ) {
      if ( line[0] == '*' )
        continue;
      std::vector<int> pdg_ids;
      std::vector<short> charges;
      double mass, mass_err_low, mass_err_high;
      double width, width_err_low, width_err_high;
      std::string part_name, part_charge_int;
      { // pdg ids
        std::istringstream ss( line.substr( 1, 32 ) );
        std::string buf;
        while ( ss >> buf )
          pdg_ids.emplace_back( std::stoi( buf ) );
      }
      { // mass + error(s)
        std::istringstream oss( line.substr( 33, 70 ) );
        oss
          >> mass >> mass_err_low >> mass_err_high;
      }
      { // width + error(s)
        std::istringstream oss( line.substr( 70, 107 ) );
        oss
          >> width >> width_err_low >> width_err_high;
      }
      { // name + charge
        std::istringstream oss( line.substr( 107 ) );
        oss
          >> part_name >> part_charge_int;
        std::istringstream oss_ch( part_charge_int );
        std::string charge_int;
        while ( std::getline( oss_ch, charge_int, ',' ) ) {
          short charge = 0;
          if ( charge_int == "-" )
            charge = -3;
          else if ( charge_int == "+" )
            charge = +3;
          else if ( charge_int == "--" )
            charge = -6;
          else if ( charge_int == "++" )
            charge = +6;
          else if ( charge_int == "0" )
            charge = 0;
          else if ( charge_int == "-1/3" )
            charge = -1;
          else if ( charge_int == "+2/3" )
            charge = +2;
          charges.emplace_back( charge );
        }
      }
      if ( pdg_ids.size() != charges.size() )
        CG_FATAL( "MCDFileParser" )
          << "Error while parsing the MCD file \"" << path << "\".";
      for ( size_t i = 0; i < pdg_ids.size(); ++i ) {
        cepgen::ParticleProperties prop{
          (cepgen::pdgid_t)pdg_ids.at( i ),
          part_name, part_name, 0,
          mass, width,
          charges.at( i ), true //FIXME
        };
        cepgen::PDG::get().define( prop );
        ++i;
      }
    }
    cepgen::PDG::get().dump();
  }
}
