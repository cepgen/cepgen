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

#include <ExRootAnalysis/ExRootTreeBranch.h>
#include <ExRootAnalysis/ExRootTreeWriter.h>
#include <TFile.h>
#include <classes/DelphesClasses.h>
#include <classes/DelphesFactory.h>
#include <modules/Delphes.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/Timer.h"

namespace cepgen {
  namespace io {
    /**
     * \brief Export handler for Delphes
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class DelphesHandler : public ExportModule {
    public:
      explicit DelphesHandler(const ParametersList&);
      ~DelphesHandler();
      static std::string description() { return "Delphes interfacing module"; }

      void initialise(const Parameters&) override;
      void setCrossSection(double cross_section, double /*err_cross_section*/) override {
        cross_section_ = cross_section;
      }
      void operator<<(const Event&) override;

    private:
      std::unique_ptr<TFile> output_;
      const std::string input_card_;
      const bool compress_;
      std::unique_ptr<Delphes> delphes_;
      //--- initialised here, but deleted elsewhere
      ExRootConfReader* conf_reader_;  // deleted at destructor
      ExRootTreeWriter* tree_writer_;  // deleted at destructor
      //--- non-owning
      DelphesFactory* factory_;
      ExRootTreeBranch* evt_branch_;
      TObjArray *out_all_parts_, *out_stab_parts_, *out_partons_;
      double cross_section_;
    };

    DelphesHandler::DelphesHandler(const ParametersList& params)
        : ExportModule(params),
          output_(new TFile(params.get<std::string>("filename", "output.delphes.root").c_str(), "recreate")),
          input_card_(params.get<std::string>("inputCard", "input.tcl")),
          compress_(params.get<bool>("compress", false)),
          delphes_(new Delphes),
          conf_reader_(new ExRootConfReader),
          tree_writer_(new ExRootTreeWriter(output_.get(), "Delphes")),
          factory_(nullptr),
          evt_branch_(nullptr),
          out_all_parts_(nullptr),
          out_stab_parts_(nullptr),
          out_partons_(nullptr),
          cross_section_(-1.) {
      try {
        conf_reader_->ReadFile(input_card_.c_str());
      } catch (const std::runtime_error& err) {
        throw CG_FATAL("DelphesHandler") << "Failed to parse the Delphes configuration card!\n\t" << err.what();
      }
      delphes_->SetTreeWriter(tree_writer_);
      delphes_->SetConfReader(conf_reader_);
    }

    DelphesHandler::~DelphesHandler() {
      delphes_->FinishTask();
      tree_writer_->Write();
    }

    void DelphesHandler::initialise(const Parameters&) {
      factory_ = delphes_->GetFactory();
      if (!factory_)
        throw CG_FATAL("DelphesHandler") << "Failed to retrieve factory object!";
      out_all_parts_ = delphes_->ExportArray("allParticles");
      out_stab_parts_ = delphes_->ExportArray("stableParticles");
      out_partons_ = delphes_->ExportArray("partons");
      evt_branch_ = tree_writer_->NewBranch("Event", LHEFEvent::Class());
      delphes_->InitTask();
    }

    void DelphesHandler::operator<<(const Event& ev) {
      delphes_->Clear();
      tree_writer_->Clear();
      //--- auxiliary event quantities
      auto evt_aux = static_cast<LHEFEvent*>(evt_branch_->NewEntry());
      evt_aux->Number = event_num_++;
      evt_aux->ProcessID = 0;
      evt_aux->Weight = ev.weight;  // events are normally unweighted in CepGen
      //evt_aux->CrossSection = cross_section_; // not yet fully supported
      evt_aux->ScalePDF = 0.;  // for the time being
      evt_aux->AlphaQED = constants::ALPHA_EM;
      evt_aux->AlphaQCD = constants::ALPHA_QCD;
      evt_aux->ReadTime = ev.time_generation;
      utils::Timer tmr;
      const auto& parts = compress_ ? ev.compress().particles() : ev.particles();
      //--- particles content
      for (const auto& part : parts) {
        auto cand = factory_->NewCandidate();
        cand->PID = part.integerPdgId();
        cand->Status = (int)part.status();
        cand->Charge = part.charge();
        //--- kinematics part
        cand->Mass = part.mass();
        const auto& mom = part.momentum();
        cand->Momentum.SetPxPyPzE(mom.px(), mom.py(), mom.pz(), mom.energy());
        // no cand->Position specified (particles produced at origin)
        //--- parentage part
        cand->M1 = part.primary() ? 0 : *part.mothers().begin();
        cand->M2 = part.mothers().size() < 2 ? 0 : *part.mothers().rbegin();
        cand->D1 = part.daughters().empty() ? -1 : *part.daughters().begin();
        cand->D2 = part.daughters().size() < 2 ? -1 : *part.daughters().rbegin();
        //--- add to the proper collection(s)
        out_all_parts_->Add(cand);
        if (cand->Status == 1)
          out_stab_parts_->Add(cand);
        else if (cand->PID <= 5 || cand->PID == 21 || cand->PID == 15)
          out_partons_->Add(cand);
      }
      delphes_->ProcessTask();
      evt_aux->ProcTime = tmr.elapsed();
      tree_writer_->Fill();
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("delphes", DelphesHandler)
