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

#include "CepGenAddOns/Herwig6Wrapper/Herwig6Interface.h"

namespace {
  extern "C" {
  void hwigin_();
  double hwuaem_(double&);
  double hwualf_(int&, double&);
  }
}  // namespace

namespace cepgen {
  namespace herwig6 {
    void initialise() {
      static bool kInitialised = false;
      if (kInitialised)
        return;
      {  // capture stdout to avoid "polluting" consumer code with unmanaged output
        int out = dup(fileno(stdout));
        freopen("/tmp/herwig.log", "w", stdout);
        hwigin_();
        dup2(out, fileno(stdout));
        close(out);
      }
      kInitialised = true;
    }
    double hwuaem(double q2) { return hwuaem_(q2); }
    double hwualf(int mode, double q2) { return hwualf_(mode, q2); }
  }  // namespace herwig6
}  // namespace cepgen
