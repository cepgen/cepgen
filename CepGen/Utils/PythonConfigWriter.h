#ifndef CepGenAddOns_BoostWrapper_PythonConfigWriter_h
#define CepGenAddOns_BoostWrapper_PythonConfigWriter_h

#include <fstream>

namespace cepgen {
  class Parameters;
  class ParametersDescription;
  namespace utils {
    class PythonConfigWriter final {
    public:
      PythonConfigWriter(const std::string&);
      ~PythonConfigWriter();

      PythonConfigWriter& operator<<(const Parameters&);
      PythonConfigWriter& operator<<(const ParametersDescription&);

    private:
      mutable std::ofstream file_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
