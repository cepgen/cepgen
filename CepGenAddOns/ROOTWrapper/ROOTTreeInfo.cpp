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
#include "CepGen/Event/Event.h"
#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

namespace ROOT {
  CepGenRun::CepGenRun() { clear(); }

  void CepGenRun::clear() {
    sqrt_s = -1.;
    xsect = errxsect = -1.;
    num_events = litigious_events = 0;
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
  }

  void CepGenEvent::fill(const cepgen::Event& ev, bool compress) {
    if (!tree_)
      throw CG_FATAL("CepGenEvent:fill") << "Trying to fill a non-existent tree!";

    clear();
    gen_time = ev.time_generation;
    tot_time = ev.time_total;
    weight = ev.weight;
    np = 0;
    const auto& parts = compress ? ev.compress().particles() : ev.particles();
    //--- loop over all particles in event
    for (const auto& part : parts) {
      const auto& mom = part.momentum();
      rapidity[np] = mom.rapidity();
      pt[np] = mom.pt();
      eta[np] = mom.eta();
      phi[np] = mom.phi();
      E[np] = part.energy();
      m[np] = part.mass();
      pdg_id[np] = part.integerPdgId();
      parent1[np] = (part.mothers().size() > 0) ? *part.mothers().begin() : -1;
      parent2[np] = (part.mothers().size() > 1) ? *part.mothers().rbegin() : -1;
      status[np] = (int)part.status();
      stable[np] = ((short)part.status() > 0);
      charge[np] = part.charge();
      role[np] = part.role();
      np++;
    }
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
    ev.time_generation = gen_time;
    ev.time_total = tot_time;
    ev.weight = weight;
    //--- first loop to populate the particles content
    for (unsigned short i = 0; i < np; ++i) {
      auto& part = ev[i];
      part.setPdgId((long)pdg_id[i]);
      part.setMomentum(cepgen::Momentum::fromPtEtaPhiE(pt[i], eta[i], phi[i], E[i]));
      part.setRole((cepgen::Particle::Role)role[i]);
      part.setStatus((cepgen::Particle::Status)status[i]);
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
