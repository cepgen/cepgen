/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include <vector>

#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main() {
  vector<cepgen::utils::Colour> colours = {cepgen::utils::Colour::none,
                                           cepgen::utils::Colour::reset,
                                           cepgen::utils::Colour::black,
                                           cepgen::utils::Colour::red,
                                           cepgen::utils::Colour::green,
                                           cepgen::utils::Colour::yellow,
                                           cepgen::utils::Colour::blue,
                                           cepgen::utils::Colour::magenta,
                                           cepgen::utils::Colour::cyan,
                                           cepgen::utils::Colour::white};
  vector<cepgen::utils::Modifier> modifiers = {cepgen::utils::Modifier::none,
                                               cepgen::utils::Modifier::reset,
                                               cepgen::utils::Modifier::bold,
                                               cepgen::utils::Modifier::dimmed,
                                               cepgen::utils::Modifier::italic,
                                               cepgen::utils::Modifier::underline,
                                               cepgen::utils::Modifier::blink,
                                               cepgen::utils::Modifier::reverse};
  CG_LOG << "Colours: " << colours;
  CG_LOG << "Modifiers: " << modifiers;
  for (const auto& col : colours) {
    CG_LOG << "<<<<< " << col << " >>>>>";
    cepgen::utils::Modifier full_mod;
    for (const auto& mod : modifiers) {
      CG_LOG << ">> " << mod << " >> " << cepgen::utils::colourise("test", col, mod);
      full_mod = full_mod | mod;
    }
    CG_LOG << "full: >>> " << cepgen::utils::colourise("test", col, full_mod);
  }
}
