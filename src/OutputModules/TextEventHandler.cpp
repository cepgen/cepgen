/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

#include <fstream>
#include <iostream>

#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Utils/Value.h"

using namespace cepgen;
using namespace std::string_literals;

/// Simple event dump module
/// \author Laurent Forthomme <laurent.forthomme@cern.ch>
/// \date Jan 2020
class TextEventHandler final : public EventExporter {
public:
  explicit TextEventHandler(const ParametersList& params)
      : EventExporter(params), save_banner_(steer<bool>("saveBanner")), print_every_(steer<int>("printEvery")) {
    if (const auto& filename = steer<std::string>("filename"); !filename.empty())
      out_ = new std::ofstream(filename);
    else
      out_ = &std::cout;
  }
  ~TextEventHandler() override {
    if (out_ != &std::cout) {
      if (auto* file_out = dynamic_cast<std::ofstream*>(out_); file_out)
        file_out->close();
      delete out_;
    }
  }

  static ParametersDescription description() {
    auto desc = EventExporter::description();
    desc.setDescription("Simple text-based event dumper");
    desc.add("saveBanner", true).setDescription("Save boilerplate in output file?");
    desc.add("printEvery", 10).setDescription("Period at which events are dumped");
    desc.add("filename", ""s).setDescription("Output filename");
    return desc;
  }

  void setCrossSection(const Value& cross_section) override {
    if (out_ != &std::cout)
      *out_ << "Total cross-section: " << cross_section << " pb.\n";
  }
  bool operator<<(const Event& ev) override {
    if (print_every_ < 0 || event_num_++ % print_every_ == 0)
      *out_ << ev << "\n";
    return true;
  }

private:
  void initialise() override {
    if (save_banner_)
      *out_ << banner("#") << "\n";
  }
  const bool save_banner_;
  const int print_every_;
  std::ostream* out_{nullptr};
};
REGISTER_EXPORTER("dump", TextEventHandler);
