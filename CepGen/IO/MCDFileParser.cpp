#include "CepGen/IO/MCDFileParser.h"
#include "CepGen/Core/Exception.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace pdg
{
  const std::unordered_map<std::string,short> MCDFileParser::m_charge_str_ = {
    { "-", -3 }, { "--", -6 },
    { "+", +3 }, { "++", +6 },
    { "0", 0 },
    { "-1/3", -1 }, { "-2/3", -2 },
    { "+1/3", +1 }, { "+2/3", +2 }
  };

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
        std::istringstream ss( line.substr( PDG_BEG, PDG_END ) );
        std::string buf;
        // split for each PDG id
        while ( ss >> buf )
          pdg_ids.emplace_back( std::stoi( buf ) );
      }
      { // mass + error(s)
        std::istringstream oss( line.substr( MASS_BEG, MASS_END ) );
        oss
          >> mass >> mass_err_low >> mass_err_high;
      }
      { // width + error(s)
        std::istringstream oss( line.substr( WIDTH_BEG, WIDTH_END ) );
        oss
          >> width >> width_err_low >> width_err_high;
      }
      { // name + charge
        std::istringstream oss( line.substr( AUX_BEG ) );
        oss
          >> part_name >> part_charge_int;
        std::istringstream oss_ch( part_charge_int );
        std::string charge_int;
        // split by ','
        while ( std::getline( oss_ch, charge_int, ',' ) ) {
          if ( m_charge_str_.count( charge_int ) == 0 )
            throw CG_FATAL( "MCDFileParser" )
              << "Failed to retrieve an integer charge "
              << "for string \"" << charge_int << "\"!";
          charges.emplace_back( m_charge_str_.at( charge_int ) );
        }
      }
      if ( pdg_ids.size() != charges.size() )
        throw CG_FATAL( "MCDFileParser" )
          << "Error while parsing the MCD file \"" << path << "\".\n\t"
          << "Invalid PDG ids / charges vectors sizes: "
          << pdg_ids.size() << " != " << charges.size() << ".";
      for ( size_t i = 0; i < pdg_ids.size(); ++i ) {
        bool is_fermion;
        short colour_ch;
        switch ( pdg_ids.at( i ) ) {
          case 1: case 2: case 3: case 4: case 5: case 6:
            colour_ch = 3;
          case 11: case 12:
          case 13: case 14:
          case 15: case 16:
            colour_ch = 0;
            is_fermion = true; break;
          default:
            colour_ch = 0;
            is_fermion = false; break;
        }
        cepgen::ParticleProperties prop{
          (cepgen::pdgid_t)pdg_ids.at( i ),
          part_name, part_name, colour_ch,
          mass, width,
          charges.at( i ), is_fermion
        };
        cepgen::PDG::get().define( prop );
        ++i;
      }
    }
    CG_INFO( "MCDFileParser" ) << "File \"" << path << "\" successfully parsed. "
      << cepgen::PDG::get().size() << " particles defined.";
  }
}
