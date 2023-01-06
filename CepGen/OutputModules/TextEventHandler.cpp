/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2023  Laurent Forthomme
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
#include <iostream>  // for cout

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"

namespace cepgen {
  namespace io {
    /// Simple event dump module
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Jan 2020
    class TextEventHandler : public ExportModule {
    public:
      explicit TextEventHandler(const ParametersList&);
      ~TextEventHandler();

      static ParametersDescription description();

      void initialise() override;
      void setCrossSection(double, double) override;
      void operator<<(const Event&) override;

    private:
      bool save_banner_;
      int print_every_;
      std::ostream* out_{nullptr};
    };

    TextEventHandler::TextEventHandler(const ParametersList& params)
        : ExportModule(params), save_banner_(steer<bool>("saveBanner")), print_every_(steer<int>("printEvery")) {
      const auto& filename = steer<std::string>("filename");
      if (!filename.empty())
        out_ = new std::ofstream(filename);
      else
        out_ = &std::cout;
    }

    TextEventHandler::~TextEventHandler() {
      if (out_ != &std::cout)
        dynamic_cast<std::ofstream*>(out_)->close();
    }

    void TextEventHandler::initialise() {
      if (save_banner_)
        *out_ << banner("#") << "\n";
    }

    void TextEventHandler::setCrossSection(double cross_section, double cross_section_err) {
      if (out_ != &std::cout)
        *out_ << "Total cross-section: " << cross_section << " +/- " << cross_section_err << " pb.\n";
    }

    void TextEventHandler::operator<<(const Event& ev) {
      if (print_every_ < 0 || event_num_++ % print_every_ == 0)
        *out_ << ev << "\n";
    }

    ParametersDescription TextEventHandler::description() {
      auto desc = ExportModule::description();
      desc.setDescription("Simple text-based event dumper");
      desc.add<bool>("saveBanner", true).setDescription("Save boilerplate in output file?");
      desc.add<int>("printEvery", 10).setDescription("Period at which events are dumped");
      desc.add<std::string>("filename", "").setDescription("Output filename");
      return desc;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("dump", TextEventHandler)
