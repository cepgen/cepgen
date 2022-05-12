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

#include <algorithm>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Hasher.h"
#include "CepGen/Utils/String.h"

extern "C" {
/// Get the particle's mass in GeV from Pythia
extern double pymass_(int&);
/// Launch the Pythia6 fragmentation
extern void pyexec_();
/// Set a parameter value to the Pythia module
extern void pygive_(const char*, int);
extern void pyckbd_();
/// List all the particles in the event in a human-readable format
extern void pylist_(int&);
/// Join two coloured particles in a colour singlet
extern void pyjoin_(int&, int&);
/// Get a particle's human-readable name from Pythia
extern void pyname_(int&, char*, int);
/// Get integer-valued event information from Pythia
extern int pyk_(int&, int&);
/// Get real-valued event information from Pythia
extern double pyp_(int&, int&);
/// Purely virtual method to call at the end of the run
void pystop_() { CG_INFO("Pythia6Hadroniser") << "End of run"; }

/// Particles content of the event
extern struct {
  /// Number of particles in the event
  int n;
  int npad;
  /// Particles' general information (status, PDG id, mother, daughter 1, daughter 2)
  int k[5][4000];
  /// Particles' kinematics, in GeV (px, py, pz, E, M)
  double p[5][4000];
  /// Primary vertex for the particles
  double v[5][4000];
} pyjets_;
}

namespace cepgen {
  namespace hadr {
    /**
     * Full interface to the Pythia 6 algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * \brief Pythia 6 hadronisation algorithm
     */
    class Pythia6Hadroniser : public Hadroniser {
    public:
      using Hadroniser::Hadroniser;

      static ParametersDescription description();

      void setRuntimeParameters(const Parameters&) override {}
      inline void readString(const char* param) override { pygive(param); }
      void init() override;
      bool run(Event& ev, double& weight, bool full) override;

      void setCrossSection(double, double) override {}

    private:
      /// Maximal number of characters to fetch for the particle's name
      static constexpr unsigned short NAME_CHR = 16;

      typedef std::unordered_map<Particle::Status, int, utils::EnumHash<Particle::Status> > ParticlesStatusMap;
      static const ParticlesStatusMap kStatusMatchMap;

      inline static double pymass(int pdgid_) { return pymass_(pdgid_); }
      inline static void pyckbd() { pyckbd_(); }
      inline static void pygive(const std::string& line) { pygive_(line.c_str(), line.length()); }
      inline static void pylist(int mlist) { pylist_(mlist); }
      inline static int pyk(int id, int qty) { return pyk_(id, qty); }
      inline static double pyp(int id, int qty) { return pyp_(id, qty); }
      inline static std::string pyname(int pdgid) {
        char out[NAME_CHR];
        std::string s;
        pyname_(pdgid, out, NAME_CHR);
        s = std::string(out, NAME_CHR);
        s.erase(remove(s.begin(), s.end(), ' '), s.end());
        return s;
      }
      /// Connect entries with colour flow information
      /// \param[in] njoin Number of particles to join in the colour flow
      /// \param[in] ijoin List of particles unique identifier to join in the colour flow
      inline static void pyjoin(std::vector<int> join) {
        int njoin = join.size();
        return pyjoin_(njoin, *join.data());
      }
      bool prepareHadronisation(Event&);
      std::pair<short, short> pickPartonsContent() const;
      size_t fillParticles(const Event&) const;
    };

    const Pythia6Hadroniser::ParticlesStatusMap Pythia6Hadroniser::kStatusMatchMap = {
        {Particle::Status::PrimordialIncoming, 21},
        {Particle::Status::FinalState, 1},
        {Particle::Status::Unfragmented, 3},
        {Particle::Status::Undecayed, 1},
        {Particle::Status::Fragmented, 11},
        {Particle::Status::Propagator, 11},
        {Particle::Status::Incoming, 11},
    };

    void Pythia6Hadroniser::init() {
      CG_WARNING("Pythia6Hadroniser") << "Branching fraction not yet implemented in this hadroniser.\n\t"
                                      << "You will have to specify manually the multiplication factor according\n\t"
                                      << "to your list of open channels.";
    }

    bool Pythia6Hadroniser::run(Event& ev, double& weight, bool full) {
      weight = 1.;

      //--- only prepare remnants for fragmentation in full (event builder) mode
      if (full && remn_fragm_)
        if (!prepareHadronisation(ev))
          return false;

      CG_DEBUG_LOOP("Pythia6Hadroniser").log([&ev](auto& dbg) {
        dbg << "Dump of the event before the hadronisation:" << ev;
      });

      //--- fill Pythia 6 common blocks
      const unsigned short str_in_evt = fillParticles(ev);

      CG_DEBUG_LOOP("Pythia6Hadroniser") << "Passed the string construction stage.\n\t "
                                         << utils::s("string object", str_in_evt, true)
                                         << " identified and constructed.";

      const int old_npart = pyjets_.n;

      //--- run the algorithm
      pyexec_();

      if (full && pyjets_.n == old_npart)
        return false;  // hadronisation failed

      //--- update the event
      for (int p = old_npart; p < pyjets_.n; ++p) {
        // filter the first particles already present in the event
        const pdgid_t pdg_id = abs(pyjets_.k[1][p]);
        ParticleProperties prop;
        if (full)
          if (!PDG::get().has(pdg_id)) {
            prop.pdgid = pdg_id;
            prop.name = pyname(pdg_id);
            prop.descr = pyname(pdg_id);
            prop.colours = pyk(p + 1, 12);  // colour factor
            prop.mass = pymass(pdg_id);
            prop.width = -1.;             //pmas( pdg_id, 2 ),
            prop.charge = pyk(p + 1, 6);  // charge
            prop.fermion = false;
            PDG::get().define(prop);
          }

        const unsigned short moth_id = pyjets_.k[2][p] - 1;
        const Particle::Role role = pyjets_.k[2][p] != 0
                                        ? ev[moth_id].role()  // child particle inherits its mother's role
                                        : Particle::Role::UnknownRole;

        auto& pa = ev.addParticle(role).get();
        pa.setId(p);
        pa.setStatus(pyjets_.k[0][p]);
        pa.setPdgId((long)pyjets_.k[1][p]);
        pa.setMomentum(Momentum(pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p]));
        pa.setMass(pyjets_.p[4][p]);
        auto& moth = ev[moth_id];
        if (role != Particle::Role::UnknownRole)
          moth.setStatus(role == Particle::Role::CentralSystem ? Particle::Status::Resonance
                                                               : Particle::Status::Fragmented);
        pa.addMother(moth);
      }
      return true;
    }

    bool Pythia6Hadroniser::prepareHadronisation(Event& ev) {
      CG_DEBUG_LOOP("Pythia6Hadroniser") << "Hadronisation preparation called.";

      for (const auto& part : ev.particles()) {
        if (part.status() != Particle::Status::Unfragmented)
          continue;
        //--- only loop over all protons to be fragmented

        const auto partons = pickPartonsContent();
        const double mx2 = part.mass2();
        const double mq = pymass(partons.first), mq2 = mq * mq;
        const double mdq = pymass(partons.second), mdq2 = mdq * mdq;

        //--- choose random direction in MX frame
        const double phi = 2. * M_PI * drand(), theta = acos(2. * drand() - 1.);  // theta angle

        //--- compute momentum of decay particles from MX
        const double px2 = 0.25 * std::pow(mx2 - mdq2 + mq2, 2) / mx2 - mq2;
        if (px2 < 0.) {
          CG_WARNING("Pythia6Hadroniser") << "Invalid remnants kinematics for " << part.role() << ".";
          return false;
        }
        const double px = std::sqrt(px2);

        //--- build 4-vectors and boost decay particles
        auto pdq = Momentum::fromPThetaPhiE(px, theta, phi, std::hypot(px, mdq));
        auto pq = -pdq;
        pq.setEnergy(std::hypot(px, mq));

        //--- singlet
        auto& quark = ev.addParticle(part.role()).get();
        quark.addMother(ev[part.id()]);
        quark.setPdgId(partons.first, +1);
        quark.setStatus(Particle::Status::Unfragmented);
        quark.setMomentum(pq.lorentzBoost(part.momentum()));

        //--- doublet
        auto& diquark = ev.addParticle(part.role()).get();
        diquark.addMother(ev[part.id()]);
        diquark.setPdgId(partons.second, +1);
        diquark.setStatus(Particle::Status::Unfragmented);
        diquark.setMomentum(pdq.lorentzBoost(part.momentum()));

        ev[part.id()].setStatus(Particle::Status::Fragmented);
      }
      return true;
    }

    size_t Pythia6Hadroniser::fillParticles(const Event& ev) const {
      //--- initialising the string fragmentation variables
      using string_t = std::vector<int>;
      std::vector<string_t> evt_strings;

      pyjets_.n = 0;  // reinitialise the event content

      for (const auto& role : ev.roles()) {  // loop on roles
        string_t evt_string;
        for (const auto& part : ev(role)) {
          const unsigned short i = part.id();
          pyjets_.p[0][i] = part.momentum().px();
          pyjets_.p[1][i] = part.momentum().py();
          pyjets_.p[2][i] = part.momentum().pz();
          pyjets_.p[3][i] = part.energy();
          pyjets_.p[4][i] = part.mass();
          try {
            pyjets_.k[0][i] = kStatusMatchMap.at(part.status());
          } catch (const std::out_of_range&) {
            ev.dump();
            throw CG_FATAL("Pythia6Hadroniser") << "Failed to retrieve a Pythia 6 particle status translation for "
                                                << "CepGen status " << (int)part.status() << "!";
          }
          pyjets_.k[1][i] = part.integerPdgId();
          const auto& moth = part.mothers();
          pyjets_.k[2][i] = moth.empty() ? 0                   // no mother
                                         : *moth.begin() + 1;  // mother
          const auto& daug = part.daughters();
          if (daug.empty())  // no daughters
            pyjets_.k[3][i] = pyjets_.k[4][i] = 0;
          else {
            pyjets_.k[3][i] = *daug.begin() + 1;   // daughter 1
            pyjets_.k[4][i] = *daug.rbegin() + 1;  // daughter 2
          }
          for (int j = 0; j < 5; ++j)
            pyjets_.v[j][i] = 0.;  // vertex position

          if (part.status() == Particle::Status::Unfragmented) {
            pyjets_.k[0][i] = 1;  // PYTHIA/JETSET workaround
            evt_string.emplace_back(part.id() + 1);
          } else if (part.status() == Particle::Status::Undecayed)
            pyjets_.k[0][i] = 2;  // intermediate resonance
          pyjets_.n++;
        }
        //--- at most one string per role
        if (!evt_string.empty())
          evt_strings.emplace_back(evt_string);
      }

      //--- loop over the strings to bind everything together
      for (const auto& evt_string : evt_strings) {
        if (evt_string.size() < 2)
          continue;

        CG_DEBUG_LOOP("Pythia6Hadroniser").log([&](auto& dbg) {
          dbg << "Joining " << utils::s("particle", evt_string.size()) << " with " << ev[evt_string[0]].role()
              << " role"
              << " in a same string";
          for (const auto& part_id : evt_string) {
            if (part_id != -1)
              dbg << utils::format("\n\t * %2d (pdgId=%4d)", part_id, pyjets_.k[1][part_id - 1]);
          }
        });
        pyjoin(evt_string);
      }
      return evt_strings.size();
    }

    std::pair<short, short> Pythia6Hadroniser::pickPartonsContent() const {
      const double ranudq = drand();
      if (ranudq < 1. / 9.)
        return {PDG::down, 2203};  // (d,uu1)
      if (ranudq < 5. / 9.)
        return {PDG::up, 2101};  // (u,ud0)
      return {PDG::up, 2103};    // (u,ud1)
    }

    ParametersDescription Pythia6Hadroniser::description() {
      auto desc = Hadroniser::description();
      desc.setDescription("Interface to the Pythia 6 string hadronisation/fragmentation algorithm");
      return desc;
    }
  }  // namespace hadr
}  // namespace cepgen

// register hadroniser
REGISTER_MODIFIER("pythia6", Pythia6Hadroniser)
