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

#ifndef CepGenAddOns_Pythia6Wrapper_Pythia6EventInterface_h
#define CepGenAddOns_Pythia6Wrapper_Pythia6EventInterface_h

#include <random>
#include <vector>

#include "CepGen/Event/Event.h"

namespace cepgen {
  namespace hadr {
    /// Interface to the Pythia 6 event content.
    class Pythia6EventInterface {
    public:
      explicit Pythia6EventInterface();

      void feedEvent(const cepgen::Event&);
      size_t numStrings() const { return evt_strings_.size(); }
      void run() const;

    private:
      /// Random number generator engine
      mutable std::default_random_engine rnd_gen_;
      std::uniform_real_distribution<double> rnd_phi_, rnd_cos_theta_;
      mutable std::uniform_real_distribution<double> rnd_qdq_;

      std::pair<short, short> pickPartonsContent() const;

      Event evt_;
      using string_t = std::vector<int>;
      std::vector<string_t> evt_strings_;
    };
  }  // namespace hadr
}  // namespace cepgen

#endif
