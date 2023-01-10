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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <ProMCBook.h>
#pragma GCC diagnostic pop

#include <cstdio>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen {
  namespace io {
    /**
     * \brief Handler for the ProMC file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class ProMCHandler : public ExportModule {
    public:
      explicit ProMCHandler(const ParametersList&);
      ~ProMCHandler();

      static ParametersDescription description();

      void initialise() override;
      void setCrossSection(double cross_section, double err) override {
        cross_section_ = cross_section, cross_section_err_ = err;
      }
      void operator<<(const Event&) override;

    private:
      static constexpr double GEV_UNIT = 1.e6;  // base unit in GEV_UNIT^-1 GeV = keV
      static constexpr double M_UNIT = 1.e3;    // base unit in M^-1 m = mm
      static int inGeV(double val) { return int(val * GEV_UNIT); }

      std::unique_ptr<ProMCBook> file_;
      const bool compress_evt_;
      const std::string log_file_path_;
      std::ofstream log_file_;
      double cross_section_{-1.}, cross_section_err_{-1.};
    };

    ProMCHandler::ProMCHandler(const ParametersList& params)
        : ExportModule(params),
          file_(new ProMCBook(steer<std::string>("filename").c_str(), "w")),
          compress_evt_(steer<bool>("compress")),
          log_file_path_(steer<std::string>("logFile")),
          log_file_(log_file_path_) {}

    ProMCHandler::~ProMCHandler() {
      ProMCStat stat;
      stat.set_cross_section_accumulated(cross_section_);
      stat.set_cross_section_error_accumulated(cross_section_err_);
      stat.set_luminosity_accumulated(event_num_ / cross_section_);
      stat.set_ntried(event_num_);
      stat.set_nselected(event_num_);
      stat.set_naccepted(event_num_);
      file_->setStatistics(stat);
      file_->close();
      //--- delete the log file once attached
      const auto num_removed_files = fs::remove_all(log_file_path_);
      CG_DEBUG("ProMCHandler") << utils::s("file", num_removed_files, true) << " removed.";
    }

    void ProMCHandler::initialise() {
      file_->setDescription(runParameters().generation().maxGen(), "Sample generated using CepGen v" + version::tag);
      log_file_ << banner() << "\n";
      ProMCHeader hdr;
      hdr.set_momentumunit(GEV_UNIT);
      hdr.set_lengthunit(M_UNIT);  // unused as for now
      for (const auto& pdg : PDG::get().particles()) {
        auto data = hdr.add_particledata();
        const auto& desc = PDG::get()(pdg);
        data->set_id((int)pdg);
        data->set_mass(desc.mass);
        data->set_name(desc.name);
        data->set_width(desc.width);
        data->set_charge(desc.charge * 1. / 3.);
      }
      hdr.set_id1(runParameters().kinematics().incomingBeams().positive().pdg);
      hdr.set_id2(runParameters().kinematics().incomingBeams().negative().pdg);
      hdr.set_pdf1(0);
      hdr.set_pdf2(0);
      hdr.set_x1(0);
      hdr.set_x2(0);
      hdr.set_ecm(runParameters().kinematics().incomingBeams().sqrtS());
      file_->setHeader(hdr);
    }

    void ProMCHandler::operator<<(const Event& ev) {
      ProMCEvent event;
      auto evt = event.mutable_event();
      evt->set_number(event_num_++);
      evt->set_process_id(0);
      evt->set_scale(ev[Particle::Role::Intermediate][0].mass());
      evt->set_alpha_qed(constants::ALPHA_EM);
      evt->set_alpha_qcd(constants::ALPHA_QCD);
      evt->set_weight(1.);

      unsigned short i = 0;
      const auto& parts = compress_evt_ ? ev.compress().particles() : ev.particles();
      for (const auto& par : parts) {
        auto part = event.mutable_particles();
        part->add_id(i++);
        part->add_pdg_id(par.integerPdgId());
        part->add_status((unsigned int)par.status());
        //--- kinematics
        part->add_px(inGeV(par.momentum().px()));
        part->add_py(inGeV(par.momentum().py()));
        part->add_pz(inGeV(par.momentum().pz()));
        part->add_energy(inGeV(par.energy()));
        part->add_mass(inGeV(par.mass()));
        part->add_barcode(0);
        //--- parentage
        const auto &daughter = par.daughters(), &moth = par.mothers();
        part->add_daughter1(daughter.empty() ? 0 : *daughter.begin() + 1);
        part->add_daughter2(daughter.size() > 1 ? *daughter.rbegin() + 1 : 0);
        part->add_mother1(moth.empty() ? 0 : *moth.begin() + 1);
        part->add_mother2(moth.size() > 1 ? *moth.rbegin() + 1 : 0);
        //--- vertex
        part->add_x(0);
        part->add_y(0);
        part->add_z(0);
        part->add_t(0);
      }
      file_->write(event);
    }

    ParametersDescription ProMCHandler::description() {
      auto desc = ExportModule::description();
      desc.setDescription("ProMC file output module");
      desc.add<std::string>("filename", "output.promc");
      desc.add<bool>("compress", false);
      desc.add<std::string>("logFile", "logfile.txt");
      return desc;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("promc", ProMCHandler)
