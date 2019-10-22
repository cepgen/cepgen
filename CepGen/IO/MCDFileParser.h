#include "CepGen/Physics/PDG.h"

namespace pdg
{
  class MCDFileParser
  {
    public:
      MCDFileParser() = default;
      static void parse( const char* path );

    private:
      static constexpr size_t PDG_BEG = 1, PDG_END = 33;
      static constexpr size_t MASS_BEG = 33, MASS_END = 70;
      static constexpr size_t WIDTH_BEG = 70, WIDTH_END = 107;
      static constexpr size_t AUX_BEG = 107;
      static const std::unordered_map<std::string,short> MAP_CHARGE_STR;
  };
}
