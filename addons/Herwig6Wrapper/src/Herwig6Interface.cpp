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

#include <unistd.h>

#include <cstdio>

#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/StreamCollector.h"
#include "CepGenHerwig6/Herwig6Interface.h"

namespace {
  extern "C" {
  void hwigin_();
  void hwsfun_(double&, double&, int&, int&, double dist[13], int&);
  double hwuaem_(double&);
  double hwualf_(int&, double&);
  }
}  // namespace

namespace cepgen::herwig6 {
  void initialise() {
    static bool kInitialised = false;
    if (kInitialised)
      return;
    std::string buf;
    {
      auto sc = utils::StreamCollector(buf);
      hwigin_();
    }
    CG_DEBUG("herwig6:initialise") << "Collected buffer at initialisation:\n" << buf;
    kInitialised = true;
  }
  double hwuaem(double q2) { return hwuaem_(q2); }
  double hwualf(int mode, double q2) { return hwualf_(mode, q2); }
  double hwsfun(double xbj, double q2, int idhad, int nset, int ibeam) {
    std::array<double, 13> dist;
    hwsfun_(xbj, q2, idhad, nset, dist.data(), ibeam);
    return dist[0];
  }
}  // namespace cepgen::herwig6
