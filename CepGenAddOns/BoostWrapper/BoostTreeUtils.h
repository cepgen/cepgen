#ifndef CepGenAddOns_BoostWrapper_BoostTreeUtils_h
#define CepGenAddOns_BoostWrapper_BoostTreeUtils_h

#include <boost/property_tree/ptree.hpp>

namespace cepgen {
  class Parameters;
  class ParametersDescription;
}  // namespace cepgen

namespace boost {
  namespace pt = property_tree;
  namespace cepgen {
    static constexpr const char* MIN_KEY = "min";
    static constexpr const char* MAX_KEY = "max";

    pt::ptree pack(const ::cepgen::Parameters&);
    pt::ptree pack(const ::cepgen::ParametersDescription&);
    pt::ptree pack(const ::cepgen::ParametersList&);
    template <typename T>
    pt::ptree pack(const std::vector<T>&);
    template <>
    pt::ptree pack(const std::vector<::cepgen::ParametersList>&);
    template <>
    pt::ptree pack(const std::vector<double>&);
    pt::ptree pack(const ::cepgen::Limits&);

    ::cepgen::ParametersList unpack(const pt::ptree&);
  }  // namespace cepgen
}  // namespace boost

#endif
