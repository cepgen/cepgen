#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include <sstream>
#include <iomanip>

namespace cepgen {
  namespace io {
    ExportModule::ExportModule(const ParametersList& params) : NamedModule(params), event_num_(0ull) {}

    ExportModule::~ExportModule() {
      CG_DEBUG("ExportModule") << "Destructor called for output module \"" << name_ << "\".";
    }

    std::string ExportModule::banner(const Parameters& params, const std::string& prep) {
      const size_t len = 45 + version::tag.size();
      std::ostringstream os;
      os << prep << "******* Sample generated with CepGen " << version::tag << " *******\n"
         << prep << " Process: " << params.processName() << " (" << params.kinematics.mode() << ")\n";
      if (params.kinematics.mode() != mode::Kinematics::ElasticElastic)
        os << prep << " Structure functions: " << params.kinematics.structureFunctions()->description() << "\n";
      if (!params.eventModifiersSequence().empty()) {
        os << prep << " " << utils::s("Event modifier", params.eventModifiersSequence().size()) << ": ";
        std::string sep;
        for (const auto& mod : params.eventModifiersSequence())
          os << sep << mod->name(), sep = ", ";
        os << "\n";
      }
      const auto& cuts = params.kinematics.cuts;
      os << prep << std::left << std::setw(len) << std::setfill('*') << "*** Incoming state "
         << "\n";
      for (const auto& cut : cuts.initial.list())
        os << prep << " " << cut.description << ": " << cut.limits << "\n";
      os << prep << std::setw(len) << std::setfill('*') << "*** Central system "
         << "\n";
      for (const auto& cut : cuts.central.list())
        os << prep << " " << cut.description << ": " << cut.limits << "\n";
      if (params.kinematics.mode() != mode::Kinematics::ElasticElastic) {
        os << prep << std::setw(len) << std::setfill('*') << "*** Remnants states "
           << "\n";
        for (const auto& cut : cuts.remnants.list())
          os << prep << " " << cut.description << ": " << cut.limits << "\n";
      }
      os << prep << std::string(45 + version::tag.size(), '*');
      return os.str();
    }
  }  // namespace io
}  // namespace cepgen
