/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
    class EventDump : public ExportModule {
    public:
      explicit EventDump(const ParametersList&);
      ~EventDump();
      static std::string description() { return "Simple text-based event dumper"; }

      void initialise(const Parameters&) override;
      void setCrossSection(double, double) override;
      void operator<<(const Event&) override;

    private:
      bool save_banner_;
      int print_every_;
      std::ostream* out_;
    };

    EventDump::EventDump(const ParametersList& params)
        : ExportModule(params),
          save_banner_(params.get<bool>("saveBanner", true)),
          print_every_(params.get<int>("printEvery", 10)),
          out_(nullptr) {
      if (params.has<std::string>("filename"))
        out_ = new std::ofstream(params.get<std::string>("filename"));
      else
        out_ = &std::cout;
    }

    EventDump::~EventDump() {
      if (out_ != &std::cout)
        dynamic_cast<std::ofstream*>(out_)->close();
    }

    void EventDump::initialise(const Parameters& params) {
      if (save_banner_)
        *out_ << banner(params, "#") << "\n";
    }

    void EventDump::setCrossSection(double cross_section, double cross_section_err) {
      if (out_ != &std::cout)
        *out_ << "Total cross-section: " << cross_section << " +/- " << cross_section_err << " pb.\n";
    }

    void EventDump::operator<<(const Event& ev) {
      if (print_every_ < 0 || event_num_++ % print_every_ == 0)
        *out_ << ev << "\n";
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("dump", EventDump)
