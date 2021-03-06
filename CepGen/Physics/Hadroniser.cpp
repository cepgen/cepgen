#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  namespace hadr {
    Hadroniser::Hadroniser(const ParametersList& plist)
        : EventModifier(plist), remn_fragm_(plist.get<bool>("remnantsFragmentation", true)) {}
  }  // namespace hadr
}  // namespace cepgen
