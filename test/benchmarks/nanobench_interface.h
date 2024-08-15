/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#ifndef CepGen_test_nanobench_interface_h
#define CepGen_test_nanobench_interface_h

#include <nanobench.h>

#include <fstream>
#include <string>
#include <vector>

#include "CepGen/Core/Exception.h"

namespace {
  inline void render_benchmark(ankerl::nanobench::Bench& bench,
                               const std::string& filename,
                               const std::vector<std::string>& outputs) {
    for (const auto& ext : outputs) {
      std::string tmpl;
      if (ext == "html")
        tmpl = ankerl::nanobench::templates::htmlBoxplot();
      else if (ext == "csv")
        tmpl = ankerl::nanobench::templates::csv();
      else if (ext == "json")
        tmpl = ankerl::nanobench::templates::json();
      else if (ext == "pyperf")
        tmpl = ankerl::nanobench::templates::pyperf();
      else
        throw CG_FATAL("main") << "Invalid output format: '" << ext << "'.";
      const auto out_filename = filename + "." + ext;
      std::ofstream out_file(out_filename);
      bench.render(tmpl, out_file);
      CG_LOG << "Successfully rendered the benchmark into '" << out_filename << "'.";
    }
  }
}  // namespace

#endif