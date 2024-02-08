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

#include "CepGen/Core/Exception.h"
#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

namespace ROOT {
  CepGenRun::CepGenRun() { clear(); }

  CepGenRun CepGenRun::load(TFile* file, const std::string& run_tree) {
    CepGenRun run;
    run.attach(file, run_tree.data());
    return run;
  }

  CepGenRun CepGenRun::load(const std::string& filename, const std::string& run_tree) {
    CepGenRun run;
    run.attach(filename.data(), run_tree.data());
    return run;
  }

  void CepGenRun::clear() {
    sqrt_s = -1.;
    xsect = errxsect = -1.;
    num_events = litigious_events = 0;
    process_name.clear();
    process_parameters.clear();
  }

  void CepGenRun::create() {
    tree_ = std::make_shared<TTree>(TREE_NAME, "a tree containing information on the previous run");
    if (!tree_)
      throw std::runtime_error("Failed to create the run TTree!");
    tree_->Branch("xsect", &xsect, "xsect/D");
    tree_->Branch("errxsect", &errxsect, "errxsect/D");
    tree_->Branch("num_events", &num_events, "num_events/i");
    tree_->Branch("litigious_events", &litigious_events, "litigious_events/i");
    tree_->Branch("sqrt_s", &sqrt_s, "sqrt_s/D");
    tree_->Branch("process_name", &process_name);
    tree_->Branch("process_parameters", &process_parameters);
  }

  void CepGenRun::fill() {
    if (!tree_)
      throw CG_FATAL("CepGenRun:fill") << "Trying to fill a non-existent tree!";
    tree_->Fill();
  }

  void CepGenRun::attach(TFile* file, const char* run_tree) {
    //--- special constructor to avoid the memory to be cleared at destruction time
    tree_ = std::shared_ptr<TTree>(dynamic_cast<TTree*>(file->Get(run_tree)), [=](TTree*) {});
    if (!tree_)
      throw std::runtime_error("Failed to attach to the run TTree!");
    tree_->SetBranchAddress("xsect", &xsect);
    tree_->SetBranchAddress("errxsect", &errxsect);
    tree_->SetBranchAddress("num_events", &num_events);
    tree_->SetBranchAddress("litigious_events", &litigious_events);
    tree_->SetBranchAddress("sqrt_s", &sqrt_s);
    auto *process_name_view = new std::string(), *process_params_view = new std::string();
    tree_->SetBranchAddress("process_name", &process_name_view);
    tree_->SetBranchAddress("process_parameters", &process_params_view);
    if (tree_->GetEntriesFast() > 1)
      std::cerr << "The run tree has more than one entry." << std::endl;
    tree_->GetEntry(0);
    process_name = *process_name_view;
    process_parameters = *process_params_view;
  }

  CepGenEvent CepGenEvent::load(TFile* file, const std::string& evt_tree) {
    CepGenEvent evt;
    evt.attach(file, evt_tree.data());
    return evt;
  }

  CepGenEvent CepGenEvent::load(const std::string& filename, const std::string& evt_tree) {
    CepGenEvent evt;
    evt.attach(filename.data(), evt_tree.data());
    return evt;
  }

  void CepGenEvent::clear() {
    gen_time = tot_time = 0.;
    np = 0;
    for (size_t i = 0; i < MAX_PART; ++i) {
      pt[i] = eta[i] = phi[i] = rapidity[i] = E[i] = m[i] = charge[i] = 0.;
      pdg_id[i] = parent1[i] = parent2[i] = stable[i] = role[i] = status[i] = 0;
    }
  }

  void CepGenEvent::create() {
    tree_ = std::make_shared<TTree>(TREE_NAME, "a tree containing information on events generated in previous run");
    if (!tree_)
      throw CG_FATAL("CepGenEvent:create") << "Failed to create the events TTree!";
    tree_->Branch("npart", &np, "npart/I");
    tree_->Branch("role", role, "role[npart]/I");
    tree_->Branch("pt", pt, "pt[npart]/D");
    tree_->Branch("eta", eta, "eta[npart]/D");
    tree_->Branch("phi", phi, "phi[npart]/D");
    tree_->Branch("rapidity", rapidity, "rapidity[npart]/D");
    tree_->Branch("E", E, "E[npart]/D");
    tree_->Branch("m", m, "m[npart]/D");
    tree_->Branch("charge", charge, "charge[npart]/D");
    tree_->Branch("pdg_id", pdg_id, "pdg_id[npart]/I");
    tree_->Branch("parent1", parent1, "parent1[npart]/I");
    tree_->Branch("parent2", parent2, "parent2[npart]/I");
    tree_->Branch("stable", stable, "stable[npart]/I");
    tree_->Branch("status", status, "status[npart]/I");
    tree_->Branch("weight", &weight, "weight/F");
    tree_->Branch("generation_time", &gen_time, "generation_time/F");
    tree_->Branch("total_time", &tot_time, "total_time/F");
    tree_->Branch("metadata", &metadata);
  }

  void CepGenEvent::fill(const cepgen::Event& ev, bool compress) {
    if (!tree_)
      throw CG_FATAL("CepGenEvent:fill") << "Trying to fill a non-existent tree!";

    clear();
    if (!ev.metadata.empty()) {
      gen_time = ev.metadata.at("time:generation");
      tot_time = ev.metadata.at("time:total");
      weight = ev.metadata.at("weight");
    }
    np = 0;
    const auto& parts = compress ? ev.compress().particles() : ev.particles();
    //--- loop over all particles in event
    for (const auto& part : parts) {
      const auto& mom = part.momentum();
      rapidity[np] = mom.rapidity();
      pt[np] = mom.pt();
      eta[np] = mom.eta();
      phi[np] = mom.phi();
      E[np] = mom.energy();
      m[np] = mom.mass();
      pdg_id[np] = part.integerPdgId();
      parent1[np] = (part.mothers().size() > 0) ? *part.mothers().begin() : -1;
      parent2[np] = (part.mothers().size() > 1) ? *part.mothers().rbegin() : -1;
      status[np] = (int)part.status();
      stable[np] = ((short)part.status() > 0);
      charge[np] = part.charge();
      role[np] = part.role();
      np++;
    }
    metadata = ev.metadata;
    tree_->Fill();
    clear();
  }

  bool CepGenEvent::next(cepgen::Event& ev) {
    if (!tree_attached_)
      try {
        attach();
      } catch (const std::runtime_error& err) {
        throw CG_FATAL("CepGenEvent:next") << "Failed to attach to the events TTree!\n" << err.what();
      }

    if (tree_->GetEntry(num_read_events_++) <= 0)
      return false;

    ev.clear();
    ev.metadata["time:generation"] = gen_time;
    ev.metadata["time:total"] = tot_time;
    ev.metadata["weight"] = weight;
    //--- first loop to populate the particles content
    for (unsigned short i = 0; i < np; ++i) {
      cepgen::Particle part;
      part.setRole((cepgen::Particle::Role)role[i]);
      part.setPdgId((long)pdg_id[i]);
      part.setStatus((cepgen::Particle::Status)status[i]);
      part.setMomentum(cepgen::Momentum::fromPtEtaPhiE(pt[i], eta[i], phi[i], E[i]));
      ev.addParticle(part);
    }
    //--- second loop to associate the parentage
    for (unsigned short i = 0; i < np; ++i) {
      auto& part = ev[i];
      if (parent1[i] > 0)
        part.addMother(ev[parent1[i]]);
      if (parent2[i] > parent1[i])
        for (unsigned short j = parent1[i] + 1; j <= parent2[i]; ++j)
          part.addMother(ev[j]);
    }
    return true;
  }
}  // namespace ROOT
