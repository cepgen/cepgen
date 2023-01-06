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

#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/Pythia8Wrapper/PythiaEventInterface.h"

namespace cepgen {
  namespace io {
    /**
     * \brief Handler for the LHE file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class LHEFPythiaHandler : public ExportModule {
    public:
      /// Class constructor
      explicit LHEFPythiaHandler(const ParametersList&);
      ~LHEFPythiaHandler();

      static ParametersDescription description();

      void initialise() override;
      /// Writer operator
      void operator<<(const Event&) override;
      void setCrossSection(double, double) override;

    private:
      std::unique_ptr<Pythia8::Pythia> pythia_;
      std::shared_ptr<Pythia8::CepGenEvent> lhaevt_;
      const bool compress_event_;
      std::string filename_;
      bool gzip_;
    };

    LHEFPythiaHandler::LHEFPythiaHandler(const ParametersList& params)
        : ExportModule(params),
          pythia_(new Pythia8::Pythia),
          lhaevt_(new Pythia8::CepGenEvent),
          compress_event_(steer<bool>("compress")),
          filename_(steer<std::string>("filename")),
          gzip_(false) {
#ifdef GZIP_BIN
      if (utils::fileExtension(filename_) == ".gz") {
        utils::replace_all(filename_, ".gz", "");
        gzip_ = true;
      }
#endif
      {
        auto file_tmp = std::ofstream(filename_);
        if (!file_tmp.is_open())
          throw CG_FATAL("LHEFPythiaHandler") << "Failed to open output filename \"" << filename_ << "\" for writing!";
      }
      lhaevt_->openLHEF(filename_);
    }

    LHEFPythiaHandler::~LHEFPythiaHandler() {
      if (lhaevt_)
        lhaevt_->closeLHEF(false);  // we do not want to rewrite the init block
#ifdef GZIP_BIN
      if (gzip_) {
        std::string cmnd(GZIP_BIN);
        cmnd += " -f " + filename_;
        if (system(cmnd.c_str()) != 0)
          CG_WARNING("LHEFPythiaHandler") << "Failed to zip the output file with command \"" << cmnd << "\".";
        CG_DEBUG("LHEFPythiaHandler") << cmnd;
      }
#endif
    }

    void LHEFPythiaHandler::initialise() {
      std::ostringstream oss_init;
      oss_init << "<!--\n" << banner() << "\n-->";
      oss_init << std::endl;  // LHEF is usually not as beautifully parsed as a standard XML...
                              // we're physicists, what do you expect?
      lhaevt_->addComments(oss_init.str());
      lhaevt_->initialise(*rt_params_);
#if PYTHIA_VERSION_INTEGER < 8300
      pythia_->setLHAupPtr(lhaevt_.get());
#else
      pythia_->setLHAupPtr(lhaevt_);
#endif
      pythia_->settings.flag("ProcessLevel:all", false);  // we do not want Pythia to interfere...
      pythia_->settings.flag("PartonLevel:all", false);   // we do not want Pythia to interfere...
      pythia_->settings.flag("HadronLevel:all", false);   // we do not want Pythia to interfere...
      pythia_->settings.mode("Beams:frameType", 5);       // LHEF event readout
      pythia_->settings.mode("Next:numberCount", 0);      // remove some of the Pythia output
      pythia_->init();
      lhaevt_->initLHEF();
    }

    void LHEFPythiaHandler::operator<<(const Event& ev) {
      lhaevt_->feedEvent(compress_event_ ? ev : ev.compress(), Pythia8::CepGenEvent::Type::centralAndFullBeamRemnants);
      pythia_->next();
      lhaevt_->eventLHEF();
    }

    void LHEFPythiaHandler::setCrossSection(double cross_section, double cross_section_err) {
      lhaevt_->setCrossSection(0, cross_section, cross_section_err);
    }

    ParametersDescription LHEFPythiaHandler::description() {
      auto desc = ExportModule::description();
      desc.setDescription("Pythia 8-based LHEF output module");
      desc.add<bool>("compress", true);
      desc.add<std::string>("filename", "output.lhe").setDescription("Output filename");
      return desc;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("lhef", LHEFPythiaHandler)
