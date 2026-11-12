/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

// CepGen includes
#include <CepGen/Core/Exception.h>
#include <CepGen/Event/Particle.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/RandomGenerator.h>
// Sherpa includes
#include <ATOOLS/Math/Vector.H>
#include <ATOOLS/Org/Exception.H>
#include <ATOOLS/Org/Return_Value.H>
#include <ATOOLS/Phys/Blob.H>
#include <ATOOLS/Phys/Blob_List.H>
#include <ATOOLS/Phys/Flavour.H>
#include <ATOOLS/Phys/Flow.H>
#include <ATOOLS/Phys/Fragmentation_Base.H>
#include <ATOOLS/Phys/Particle.H>
#include <SHERPA/Single_Events/Event_Phase_Handler.H>

#include "CepGenSherpa/DiffractiveState.h"
#include "CepGenSherpa/SherpaFramework.h"

namespace cepgen::sherpa {
  namespace {
    constexpr size_t kMaxSteps = 1000;  // maximal number of Sherpa event phases calls for a single system

    /// Build an on-shell parton of the colour-singlet system, in its rest frame
    ATOOLS::Particle* buildParton(long pdg_id, const Momentum& p3, double mass, int colour, bool spectator) {
      const ATOOLS::Flavour flavour(static_cast<kf_code>(std::labs(pdg_id)), pdg_id < 0);
      const auto energy = std::hypot(p3.p(), mass);
      auto* parton = new ATOOLS::Particle(
          -1, flavour, ATOOLS::Vec4D(energy, p3.px(), p3.py(), p3.pz()), spectator ? 'B' : 'F');  // garbage-collected
      parton->SetNumber();
      parton->SetStatus(ATOOLS::part_status::fragmented);
      parton->SetFinalMass(mass);
      // a (anti)quark-(anti)diquark pair is a colour singlet: one colour (triplet) and one anticolour (antitriplet) line
      if (flavour.StrongCharge() > 0)
        parton->SetFlow(1, colour);
      else
        parton->SetFlow(2, colour);
      if (spectator)
        parton->SetBeam(0);  // flagged as a beam remnant, for the specific treatment of the fragmentation algorithm
      return parton;
    }

    /// Make a Sherpa hadron known to the CepGen particles table
    void registerParticle(const ATOOLS::Flavour& flavour) {
      const auto pdg_id = static_cast<pdgid_t>(flavour.Kfcode());
      if (PDG::get().has(pdg_id))
        return;
      const ATOOLS::Flavour particle(flavour.Kfcode());  // particle (not antiparticle) properties
      ParticleProperties properties;
      properties.pdgid = pdg_id;
      properties.name = particle.IDName();
      properties.human_name = particle.IDName();
      properties.colours = 1;
      properties.mass = particle.Mass(true);
      properties.width = particle.Width();
      if (const auto charge = particle.IntCharge(); charge != 0)  // in units of e/3 for both Sherpa and CepGen
        properties.charges = {charge, -charge};
      properties.fermion = particle.IsFermion();
      PDG::get().define(properties);
    }

    /// Run the secondary event phases until none of them has anything left to treat, as in the Sherpa event handler
    bool runSecondaryPhases(ATOOLS::Blob_List& blobs, const Framework& framework) {
      const auto& phases = framework.secondaryPhases();
      size_t num_steps = 0;
      for (auto phase = phases.begin(); phase != phases.end();) {
        if (++num_steps > kMaxSteps) {
          CG_WARNING("sherpa:hadronise") << "Maximal number of Sherpa event phases calls (" << kMaxSteps
                                         << ") reached.";
          return false;
        }
        switch (const auto result = (*phase)->Treat(&blobs); result) {
          case ATOOLS::Return_Value::Success:  // something changed, all phases have to be checked again
            phase = phases.begin();
            break;
          case ATOOLS::Return_Value::Nothing:
            ++phase;
            break;
          case ATOOLS::Return_Value::Retry_Phase:
            break;
          default:
            CG_DEBUG("sherpa:hadronise") << "Sherpa phase '" << (*phase)->Name() << "' returned " << result << ".";
            return false;
        }
      }
      return true;
    }

    /// Collect all colourless particles, following the decay chains from the system fragmentation
    void collect(ATOOLS::Blob* blob,
                 int mother,
                 hadr::Hadroniser::Fragments& fragments,
                 std::unordered_set<const ATOOLS::Blob*>& visited_blobs) {
      if (!blob || !visited_blobs.insert(blob).second)  // multi-partons blobs are reached once per incoming parton
        return;
      for (int i = 0; i < blob->NOutP(); ++i) {
        auto* particle = blob->OutParticle(i);
        auto* decay = particle->DecayBlob();
        auto parent = mother;
        if (particle->Flav().StrongCharge() == 0) {  // intermediate partons are not exported, only their products
          registerParticle(particle->Flav());
          hadr::Hadroniser::Fragment fragment;
          fragment.index = static_cast<int>(fragments.size());
          fragment.mother = mother;
          fragment.pdg_id = static_cast<long>(particle->Flav());
          fragment.stable = !decay || decay->NOutP() == 0;
          const auto& p = particle->Momentum();
          fragment.momentum = Momentum::fromPxPyPzE(p[1], p[2], p[3], p[0]);
          fragments.emplace_back(fragment);
          parent = fragment.index;
        }
        collect(decay, parent, fragments, visited_blobs);
      }
    }
  }  // namespace

  DiffractiveState::DiffractiveState(const Particle& beam,
                                     const Momentum& system,
                                     const std::function<double(long, bool)>& parton_mass,
                                     utils::RandomGenerator& rng,
                                     double primordial_kt)
      : hadr::Hadroniser::DiffractiveState(system), system_pdg_id_(beam.integerPdgId()) {
    if (std::labs(system_pdg_id_) != PDG::proton)
      throw CG_FATAL("sherpa:DiffractiveState")
          << "Beam particle with PDG id " << system_pdg_id_
          << " is not supported. Only protons and antiprotons can be dissociated.";

    // valence flavour content, weighted by the quark electric charge squared and its multiplicity in the proton
    const long sign = system_pdg_id_ < 0 ? -1 : +1;
    if (rng.uniform() < 1. / 9.) {  // valence down quark
      valquark_pdg_id_ = sign * PDG::down;
      diquark_pdg_id_ = sign * PDG::uu1;
    } else {  // valence up quark
      valquark_pdg_id_ = sign * PDG::up;
      diquark_pdg_id_ = sign * PDG::ud0;
    }
    for (const bool constituent : {true, false}) {  // fall back to perturbative masses for the lightest systems
      constituent_masses_ = constituent;
      valquark_mass_ = parton_mass(valquark_pdg_id_, constituent);
      diquark_mass_ = parton_mass(diquark_pdg_id_, constituent);
      if (valquark_mass_ + diquark_mass_ < system_mass_)
        break;
    }

    const auto mass2 = system_mass_ * system_mass_;
    const auto sum = valquark_mass_ + diquark_mass_, diff = valquark_mass_ - diquark_mass_;
    const auto lambda = (mass2 - sum * sum) * (mass2 - diff * diff);
    if (system_mass_ <= sum || lambda <= 0.)
      throw CG_FATAL("sherpa:DiffractiveState")
          << "Dissociated system of mass " << system_mass_ << " GeV is too light for its constituents.";
    const auto p_star = 0.5 * std::sqrt(lambda) / system_mass_;

    // the spectator follows the beam direction, as seen from the rest frame of the dissociated system
    const auto beam_in_system =
        Momentum(beam.momentum())
            .lorentzBoost(Momentum::fromPxPyPzE(-system.px(), -system.py(), -system.pz(), system.energy()));
    const auto axis =
        beam_in_system.p() > 0. ? beam_in_system.unit() : Momentum::randomThetaPhi(rng, {-M_PI, M_PI}, {0., 2. * M_PI});

    auto direction = axis;
    if (primordial_kt > 0.) {
      const auto px = rng.normal(0., primordial_kt * M_SQRT1_2), py = rng.normal(0., primordial_kt * M_SQRT1_2);
      if (const auto pt = std::hypot(px, py); pt > 0.) {
        const auto sin_theta = std::min(pt / p_star, 1.), cos_theta = std::sqrt(1. - sin_theta * sin_theta);
        const auto helper =
            std::fabs(axis.px()) < 0.9 ? Momentum::fromPxPyPzE(1., 0., 0., 1.) : Momentum::fromPxPyPzE(0., 1., 0., 1.);
        const auto e1 = axis.crossProduct(helper).unit(), e2 = axis.crossProduct(e1);
        direction = cos_theta * axis + sin_theta * (px / pt * e1 + py / pt * e2);
      }
    }
    diquark_momentum_ = p_star * direction;
  }

  bool DiffractiveState::hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const {
    {  // start by preparing the environment after any previous run
      fragments.clear();
      ATOOLS::Blob::Reset();
      ATOOLS::Particle::Reset();
      ATOOLS::Flow::ResetCounter();
    }
    ATOOLS::Blob_List blobs(/*destruct=*/true);  // owns the blobs, and therefore all particles attached to them
                                                 // -> garbage-collection

    // the system is hadronised at rest, as a colour singlet ready for fragmentation; the Sherpa algorithms expect the
    // partons to be attached to a full chain of upstream blobs:
    //   [bunch: beam -> system] -> [beam: system -> partons] -> [fragmentation: partons -> hadrons]
    const ATOOLS::Flavour system_flavour(static_cast<kf_code>(std::labs(system_pdg_id_)), system_pdg_id_ < 0);
    const ATOOLS::Vec4D system_momentum(system_mass_, 0., 0., 0.);
    auto* incoming_system = new ATOOLS::Particle(-1, system_flavour, system_momentum, 'B');
    auto* system = new ATOOLS::Particle(-1, system_flavour, system_momentum, 'B');  // garbage-collected
    for (auto* particle : {incoming_system, system}) {
      particle->SetNumber();
      particle->SetFinalMass(system_mass_);
      particle->SetStatus(ATOOLS::part_status::decayed);
    }
    auto* bunch_blob = new ATOOLS::Blob();
    bunch_blob->SetType(ATOOLS::btp::Bunch);
    bunch_blob->SetTypeSpec("CepGen_dissociated_beam");
    bunch_blob->SetId();
    bunch_blob->AddToInParticles(incoming_system);
    bunch_blob->AddToOutParticles(system);
    blobs.push_back(bunch_blob);  // garbage-collected

    auto* remnants_blob = new ATOOLS::Blob();
    remnants_blob->SetType(ATOOLS::btp::Beam);
    remnants_blob->SetTypeSpec("CepGen_dissociated_system");
    remnants_blob->SetId();
    remnants_blob->AddToInParticles(system);
    blobs.push_back(remnants_blob);  // garbage-collected

    auto* blob = new ATOOLS::Blob();
    blob->SetType(ATOOLS::btp::Fragmentation);
    blob->SetTypeSpec("CepGen_diffractive_system");
    blob->SetStatus(ATOOLS::blob_status::needs_hadronization);
    blob->SetId();
    const auto colour = static_cast<int>(ATOOLS::Flow::Counter());
    for (auto* parton : {buildParton(valquark_pdg_id_, -diquark_momentum_, valquark_mass_, colour, false),
                         buildParton(diquark_pdg_id_, diquark_momentum_, diquark_mass_, colour, true)}) {
      remnants_blob->AddToOutParticles(parton);
      blob->AddToInParticles(parton);
    }
    blobs.push_back(blob);  // garbage-collected

    bool success = false;
    try {  // run fragmentation
      // system blob fed to fragmentation algorithm, without:
      //  - singlets sorting (only one colour singlet),
      //  - colour reconnections (idem),
      //  - momentum stretching (allows systems lighter than sum(m_const))
      if (const auto result = Framework::get().fragmentation().Hadronize(&blobs);
          result != ATOOLS::Return_Value::Success && result != ATOOLS::Return_Value::Nothing) {
        CG_DEBUG("sherpa:hadronise") << "Fragmentation failed with return code " << result << ".";
      } else if (!blob->Has(ATOOLS::blob_status::needs_hadronization) && blob->NOutP() > 0)
        success = runSecondaryPhases(blobs, Framework::get());
    } catch (const ATOOLS::Return_Value::code& code) {
      CG_DEBUG("sherpa:hadronise") << "Sherpa requested a new trial (return code " << code << ").";
    } catch (const ATOOLS::Exception& exception) {
      CG_WARNING("sherpa:hadronise").log([&exception](auto& log) { log << "Sherpa exception caught: " << exception; });
    }
    for (auto* phase : Framework::get().secondaryPhases())
      phase->CleanUp();
    if (!success)
      return false;

    std::unordered_set<const ATOOLS::Blob*> visited_blobs;
    collect(blob, -1, fragments, visited_blobs);
    if (!fragments.empty())
      return momentumConserved(fragments, tolerance);
    return false;
  }
}  // namespace cepgen::sherpa
