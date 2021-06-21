#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  EventModifier::EventModifier(const ParametersList& plist)
      : NamedModule(plist),
        seed_(plist.get<int>("seed", -1)),
        max_trials_(plist.get<int>("maxTrials", 1)),
        rt_params_(nullptr) {
    CG_DEBUG("EventModifier:init") << "\"" << name_ << "\"-type event modifier built with:\n\t"
                                   << "* seed = " << seed_ << "\n\t"
                                   << "* maximum trials: " << max_trials_;
  }

  void EventModifier::readStrings(const std::vector<std::string>& params) {
    if (params.empty())
      return;
    std::ostringstream os;
    for (const auto& p : params) {
      readString(p);
      os << "\n\t  '" << p << "'";
    }
    CG_DEBUG("EventModifier:configure") << "Feeding \"" << name_ << "\" event modifier algorithm with:" << os.str();
  }
}  // namespace cepgen
