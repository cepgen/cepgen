/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Value.h"
#include "CepGen/Version.h"

using namespace cepgen;

/// Handler for the ProMC file output
/// \author Laurent Forthomme <laurent.forthomme@cern.ch>
/// \date Jul 2019
class ProMCHandler final : public EventExporter {
public:
  explicit ProMCHandler(const ParametersList& params)
      : EventExporter(params),
        compress_evt_(steer<bool>("compress")),
        log_file_path_(steer<std::string>("logFile")),
        log_file_(log_file_path_) {}

  ~ProMCHandler() override {
    if (file_) {
      ProMCStat stat;
      stat.set_cross_section_accumulated(cross_section_);
      stat.set_cross_section_error_accumulated(cross_section_.uncertainty());
      stat.set_luminosity_accumulated(event_num_ / cross_section_);
      stat.set_ntried(event_num_);
      stat.set_nselected(event_num_);
      stat.set_naccepted(event_num_);
      file_->setStatistics(stat);
      file_->close();
    }
    const auto num_removed_files = fs::remove_all(log_file_path_);  // delete the log file once attached
    CG_DEBUG("ProMCHandler") << utils::s("file", num_removed_files, true) << " removed.";
  }

  static ParametersDescription description() {
    auto desc = EventExporter::description();
    desc.setDescription("ProMC file output module");
    desc.add<std::string>("filename", "output.promc");
    desc.add<bool>("compress", false);
    desc.add<std::string>("logFile", "logfile.txt");
    return desc;
  }

  void setCrossSection(const Value& cross_section) override { cross_section_ = cross_section; }
  bool operator<<(const Event& event) override {
    if (!file_)
      return false;
    ProMCEvent promc_event;
    auto evt = promc_event.mutable_event();
    evt->set_number(event_num_++);
    evt->set_process_id(0);
    evt->set_scale(event.oneWithRole(Particle::Role::Intermediate).momentum().mass());
    evt->set_alpha_qed(event.metadata("alphaEM"));
    evt->set_alpha_qcd(event.metadata("alphaS"));
    evt->set_weight(1.);

    unsigned short i = 0;
    const auto& parts = compress_evt_ ? event.compress().particles() : event.particles();
    for (const auto& par : parts) {
      auto part = promc_event.mutable_particles();
      part->add_id(i++);
      part->add_pdg_id(par.integerPdgId());
      part->add_status(static_cast<unsigned int>(par.status()));
      //--- kinematics
      part->add_px(inGeV(par.momentum().px()));
      part->add_py(inGeV(par.momentum().py()));
      part->add_pz(inGeV(par.momentum().pz()));
      part->add_energy(inGeV(par.momentum().energy()));
      part->add_mass(inGeV(par.momentum().mass()));
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
    return file_->write(promc_event);
  }

private:
  void initialise() override {
    file_.reset(new ProMCBook(steer<std::string>("filename").c_str(), "w"));
    file_->setDescription(runParameters().generation().maxGen(), "Sample generated using CepGen v" + version::tag);
    log_file_ << banner() << "\n";
    ProMCHeader hdr;
    hdr.set_momentumunit(GEV_UNIT);
    hdr.set_lengthunit(M_UNIT);  // unused as for now
    for (const auto& pdg : PDG::get().particles()) {
      auto data = hdr.add_particledata();
      const auto& desc = PDG::get()(pdg);
      data->set_id(static_cast<int>(pdg));
      data->set_mass(desc.mass);
      data->set_name(desc.name);
      data->set_width(desc.width);
      data->set_charge(desc.integerCharge() * 1. / 3.);
    }
    hdr.set_id1(runParameters().kinematics().incomingBeams().positive().integerPdgId());
    hdr.set_id2(runParameters().kinematics().incomingBeams().negative().integerPdgId());
    hdr.set_pdf1(0);
    hdr.set_pdf2(0);
    hdr.set_x1(0);
    hdr.set_x2(0);
    hdr.set_ecm(runParameters().kinematics().incomingBeams().sqrtS());
    file_->setHeader(hdr);
  }

  static constexpr double GEV_UNIT = 1.e6;  // base unit in GEV_UNIT^-1 GeV = keV
  static constexpr double M_UNIT = 1.e3;    // base unit in M^-1 m = mm
  static int inGeV(double val) { return static_cast<int>(val * GEV_UNIT); }

  std::unique_ptr<ProMCBook> file_;
  const bool compress_evt_;
  const std::string log_file_path_;
  std::ofstream log_file_;
  Value cross_section_{0., 1.};
};
REGISTER_EXPORTER("promc", ProMCHandler);
