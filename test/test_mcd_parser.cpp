#include "CepGen/IO/MCDFileParser.h"

int main()
{
  pdg::MCDFileParser::parse( "External/mass_width_2019.mcd" );
  cepgen::PDG::get().dump();
  return 0;
}
