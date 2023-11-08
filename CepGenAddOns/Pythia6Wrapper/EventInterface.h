/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#ifndef CepGenAddOns_Pythia6Wrapper_EventInterface_h
#define CepGenAddOns_Pythia6Wrapper_EventInterface_h

#include <random>
#include <vector>

namespace cepgen {
  class Event;
}
namespace pythia6 {
  /// Interface to the Pythia 6 event content.
  class EventInterface {
  public:
    explicit EventInterface(cepgen::Event&);

    void prepareHadronisation();
    size_t numStrings() const { return evt_strings_.size(); }
    void run();

  private:
    void fillEventBlock();

    cepgen::Event& evt_;  // NOT owning

    std::default_random_engine rnd_gen_;  ///< Random number generator engine
    std::uniform_real_distribution<double> rnd_phi_, rnd_cos_theta_, rnd_qdq_;

    std::pair<short, short> pickPartonsContent();

    using string_t = std::vector<int>;
    std::vector<string_t> evt_strings_;
  };
}  // namespace pythia6

#endif
