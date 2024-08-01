/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#ifndef CepGen_Utils_Collections_h
#define CepGen_Utils_Collections_h

#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cepgen::utils {
  /// Retrieve all keys from a map
  template <typename K, typename T>
  inline std::vector<K> keys(const std::map<K, T>& coll) {
    std::vector<K> keys;
    std::transform(coll.begin(), coll.end(), std::back_inserter(keys), [](const auto& val) { return val.first; });
    return keys;
  }
  /// Retrieve all keys from an unordered map
  template <typename K, typename T>
  inline std::vector<K> keys(const std::unordered_map<K, T>& coll) {
    std::vector<K> keys;
    std::transform(coll.begin(), coll.end(), std::back_inserter(keys), [](const auto& val) { return val.first; });
    return keys;
  }
  /// Check if a vector contains an item
  template <typename T>
  inline bool contains(const std::vector<T>& coll, const T& item) {
    return std::find(coll.begin(), coll.end(), item) != coll.end();
  }
  /// Check if a set contains an item
  template <typename T>
  inline bool contains(const std::set<T>& coll, const T& item) {
    return std::find(coll.begin(), coll.end(), item) != coll.end();
  }
  /// Check if an unordered map contains an item
  template <typename K, typename T>
  inline bool contains(const std::unordered_map<K, T>& coll, const T& item) {
    return std::find_if(coll.begin(), coll.end(), [&item](const auto& kv) { return kv.second == item; }) != coll.end();
  }
  /// Remove duplicates and sort a collection
  template <typename T>
  inline void normalise(std::vector<T>& coll) {
    std::unordered_set<T> set;
    for (const auto& it : coll)
      set.insert(it);
    coll.assign(set.begin(), set.end());
    std::sort(coll.begin(), coll.end());
  }
  /// Check if all elements of a collection are uniform
  template <typename T>
  inline bool uniform(const std::vector<T>& coll) {
    return coll.size() > 1 ? coll == std::vector<T>(coll.size(), coll.at(0)) : true;
  }
}  // namespace cepgen::utils

#endif
