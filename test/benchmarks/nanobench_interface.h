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
