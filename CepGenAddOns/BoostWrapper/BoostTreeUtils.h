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

  namespace boost {
    namespace pt = ::boost::property_tree;
    static constexpr const char* MIN_KEY = "min";
    static constexpr const char* MAX_KEY = "max";
    static constexpr const char* DAUGH_KEY = "DAUGHTER";

    pt::ptree pack(const Parameters&);
    pt::ptree pack(const ParametersDescription&);
    pt::ptree pack(const ParametersList&);
    template <typename T>
    pt::ptree pack(const std::vector<T>&);
    template <>
    pt::ptree pack(const std::vector<ParametersList>&);
    template <>
    pt::ptree pack(const std::vector<double>&);
    pt::ptree pack(const Limits&);

    ParametersList unpack(const pt::ptree&);
  }  // namespace boost
}  // namespace cepgen

#endif
