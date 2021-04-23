#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Core/Exception.h"

namespace ROOT {
  void CepGenEvent::fill(const cepgen::Event& ev, bool compress) {
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
    fill();
  }

  bool CepGenEvent::next(cepgen::Event& ev) {
    if (!tree_attached_)
      try {
        attach();
      } catch (const std::runtime_error& err) {
        throw CG_FATAL("ROOT::CepGenEvent") << "Failed to attach to the events TTree!\n" << err.what();
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
