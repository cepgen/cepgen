#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/PythonConfigWriter.h"

namespace cepgen {
  namespace utils {
    PythonConfigWriter::PythonConfigWriter(const std::string& filename) : file_(filename) {
      file_ << "import Config.Core as cepgen\n\n";
    }

    PythonConfigWriter::~PythonConfigWriter() { file_.close(); }

    PythonConfigWriter& PythonConfigWriter::operator<<(const Parameters&) { return *this; }

    PythonConfigWriter& PythonConfigWriter::operator<<(const ParametersDescription& pdesc) {
      CG_DEBUG("PythonConfigWriter") << "Adding a parameters description object:\n" << pdesc;
      std::function<std::string(const ParametersDescription&, const std::string&, size_t)> write =
          [&](const ParametersDescription& pdesc, const std::string& key, size_t offset) -> std::string {
        // write a generic parameters description object
        std::stringstream os;
        const std::string off(offset * 4, ' ');
        os << off;
        if (!key.empty())
          os << key << " = ";

        const auto& params = pdesc.parameters();
        switch (pdesc.type()) {
          case ParametersDescription::Type::Module:
            os << "cepgen.Module(" << params.getString(ParametersList::MODULE_NAME) << ",";
            break;
          case ParametersDescription::Type::Parameters:
            os << "cepgen.Parameters(";
            break;
          case ParametersDescription::Type::Value:
            break;
        }
        std::string sep = "";
        for (const auto& key : params.keys(false)) {
          os << sep << "\n";
          const auto& daugh = pdesc.get(key);
          switch (daugh.type()) {
            case ParametersDescription::Type::Module:
            case ParametersDescription::Type::Parameters:
              os << write(pdesc.get(key), key, offset + 1);
              break;
            case ParametersDescription::Type::Value:
              os << off << std::string(4, ' ') << key << " = " << params.getString(key, true);
              break;
          }
          sep = ",";
        }
        switch (pdesc.type()) {
          case ParametersDescription::Type::Module:
          case ParametersDescription::Type::Parameters:
            os << "\n" << off;
            break;
          case ParametersDescription::Type::Value:
            break;
        }
        os << ")";
        return os.str();
      };

      file_ << write(pdesc, "", 0) << "\n";

      return *this;
    }
  }  // namespace utils
}  // namespace cepgen
