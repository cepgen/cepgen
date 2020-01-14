#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"

int main()
{
  pdg::MCDFileParser::parse( "External/mass_width_2019.mcd" );
  cepgen::PDG::get().dump();
  return 0;
}
