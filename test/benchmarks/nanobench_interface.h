#include <nanobench.h>

#include <fstream>
#include <string>
#include <vector>

#include "CepGen/Core/Exception.h"

namespace {
  inline void render_benchmark(ankerl::nanobench::Bench& bench, const std::vector<std::string>& outputs) {
    for (const auto& output : outputs) {
      std::string tmpl;
      if (output == "html")
        tmpl = ankerl::nanobench::templates::htmlBoxplot();
      else if (output == "csv")
        tmpl = ankerl::nanobench::templates::csv();
      else if (output == "json")
        tmpl = ankerl::nanobench::templates::json();
      else if (output == "pyperf")
        tmpl = ankerl::nanobench::templates::pyperf();
      else
        throw CG_FATAL("main") << "Invalid output format: '" << output << "'.";
      std::ofstream out_file("benchmark." + output);
      bench.render(tmpl, out_file);
    }
  }
}  // namespace
