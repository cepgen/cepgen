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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportHandler.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Parameters.h"
#include "CepGenPGS/PGSInterface.h"

namespace cepgen::pgs {
  /// PGS export handler
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Aug 2019
  class ExportModule : public cepgen::ExportModule {
  public:
    explicit ExportModule(const ParametersList& params)
        : cepgen::ExportModule(params),
          compress_(steer<bool>("compress")),
          do_trig_(steer<bool>("simulateTrigger")),
          do_reco_(steer<bool>("simulateReco")) {
      steer<std::string>("inputCard").copy(pgsevt_.pgs_param_file, sizeof(pgsevt_.pgs_param_file));
      steer<std::string>("outputFile").copy(pgsevt_.pgs_output_file, sizeof(pgsevt_.pgs_output_file));
      //steer<std::string>("outputLog").copy(pgsevt_.pgs_log_file, sizeof(pgsevt_.pgs_log_file));  // unused
      pgsevt_.pgs_log_unit = 6;  // "stdout"
      pgsevt_.numarg = 0;
      std::string("events").copy(pgsevt_.evtlum.data(), pgsevt_.size());
      std::string("USER").copy(pgsevt_.optpgs.data(), pgsevt_.optpgs.size());
      if (do_trig_)  // calculate mask for printout operation
        mask_ += pgs::printmask::trg_obj;
      if (do_reco_) {
        //mask_ += pgs::printmask::calo_sum;
        //mask_ += pgs::printmask::calo_clus;
        mask_ += pgs::printmask::off_obj;
      }
    }
    ~ExportModule() { pgs_write_event_("end", 3); }

    ParametersDescription description() {
      auto desc = cepgen::ExportModule::description();
      desc.add<bool>("compress", true);
      desc.add<bool>("simulateTrigger", true);
      desc.add<bool>("simulateReco", false);
      desc.add<std::string>("inputCard", "lhc.par");
      desc.add<std::string>("outputFile", "pgs.out");
      desc.add<std::string>("outputLog", "pgs_out.log");
      return desc;
    }

    void initialise(const Parameters& params) override {
      pgsevt_.nevpgs = params.generation().maxgen;
      pgsevt_.nprpgs = params.generation().gen_print_every;  // in PGS, only used for fragmentation algos
                                                             // reusing it here to save a private variable
      pgs_initialize_();
      pgs_write_event_("begin", 5);
      CG_INFO("pgs:ExportModule") << "PGS initialised with input parameters card at:\n\n\t"
                                  << "  \"" << pgsevt_.pgs_param_file << "\".\n\n\t"
                                  << "Trigger emulation: " << do_trig_ << "\n\t"
                                  << "Detector reconstruction: " << do_reco_ << ".";
    }
    void operator<<(const Event& ev) override {
      hepevt_.nevhep = ++event_num_;
      const auto& parts = compress_ ? ev.compressed().particles() : ev.particles();
      hepevt_.nhep = 0;
      //--- particles content
      for (const auto& part : parts) {
        hepevt_.isthep[hepevt_.nhep] = (int)part.status();
        hepevt_.idhep[hepevt_.nhep] = part.integerPdgId();
        hepevt_.jmohep[hepevt_.nhep][0] = part.primary() ? 0 : *part.mothers().begin() + 1;
        hepevt_.jmohep[hepevt_.nhep][1] = part.mothers().size() < 2 ? 0 : *part.mothers().rbegin() + 1;
        hepevt_.jdahep[hepevt_.nhep][0] = part.daughters().empty() ? 0 : *part.daughters().begin() + 1;
        hepevt_.jdahep[hepevt_.nhep][1] = part.daughters().size() < 2 ? 0 : *part.daughters().rbegin() + 1;
        const auto& mom = part.momentum().pVector();
        std::copy(mom.begin(), mom.end(), hepevt_.phep[hepevt_.nhep]);
        for (unsigned short i = 0; i < 4; ++i)
          hepevt_.vhep[hepevt_.nhep][i] = 0.;
        hepevt_.nhep++;
      }
      pgs_next_event_();

      if (do_trig_)  // perform trigger simulation
        pgs_trigger_();
      if (do_reco_)  // perform reconstruction; results stored in common block
        pgs_recon_();
      //if (hepevt_.nevhep % pgsevt_.nprpgs == 0) {
      double hepcut = 0., calcut = 0.;
      pgs_dump_event_(mask_, hepcut, calcut);
      CG_INFO("") << pgsrec_.numobj;
      //}
      pgs_write_event_("event", 5);
    }

  private:
    const bool compress_;
    const bool do_trig_, do_reco_;
    int mask_{pgs::printmask::hepevt};
  };
}  // namespace cepgen::pgs
using PGSExportModule = cepgen::pgs::ExportModule;
REGISTER_IO_MODULE("pgs", PGSExportModule)
