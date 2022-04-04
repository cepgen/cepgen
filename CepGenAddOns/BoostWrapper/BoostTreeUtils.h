/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
