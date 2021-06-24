#include "CepGenAddOns/Pythia8Wrapper/PythiaEventInterface.h"

#include "CepGen/Core/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Filesystem.h"

#include "Pythia8/Pythia.h"

#include <sstream>

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
      static std::string description() { return "Pythia 8-based LHEF output module"; }

      void initialise(const Parameters&) override;
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
          compress_event_(params.get<bool>("compress", true)),
          filename_(params.get<std::string>("filename", "output.lhe")),
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
    }  // namespace io

    LHEFPythiaHandler::~LHEFPythiaHandler() {
      if (lhaevt_)
        lhaevt_->closeLHEF(false);  // we do not want to rewrite the init block
#ifdef GZIP_BIN
      if (gzip_) {
        std::string cmnd(GZIP_BIN);
        cmnd += " -f " + filename_;
        system(cmnd.c_str());
        CG_INFO("") << cmnd;
      }
#endif
    }

    void LHEFPythiaHandler::initialise(const Parameters& params) {
      std::ostringstream oss_init;
      oss_init << "<!--\n" << banner(params) << "\n-->";
      oss_init << std::endl;  // LHEF is usually not as beautifully parsed as a standard XML...
                              // we're physicists, what do you expect?
      lhaevt_->addComments(oss_init.str());
      lhaevt_->initialise(params);
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
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("lhef", LHEFPythiaHandler)
