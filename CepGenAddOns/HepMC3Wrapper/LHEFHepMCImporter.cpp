/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include <HepMC3/LHEF.h>
#include <HepMC3/Version.h>

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/HepMC3Wrapper/HepMC3EventInterface.h"

namespace cepgen {
  /// HepMC3 handler for the LHEF file import
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Feb 2024
  class LHEFHepMCImporter final : public EventImporter {
  public:
    /// Class constructor
    explicit LHEFHepMCImporter(const ParametersList& params) : EventImporter(params) {
      if (const auto filename = steer<std::string>("filename"); !filename.empty()) {
        try {
          reader_.reset(new LHEF::Reader(filename));
        } catch (std::runtime_error& err) {
          throw CG_FATAL("LHEFHepMCImporter") << "Failed to load the LHEF file. Error:\n" << err.what();
        }
        CG_INFO("LHEFHepMCImporter") << "Interfacing module initialised "
                                     << "for HepMC version " << HEPMC3_VERSION << " and LHEF file '" << filename
                                     << "' with version " << reader_->version << ".";
      } else
        throw CG_FATAL("LHEFHepMCImporter") << "Failed to retrieve the file name from module builder attributes.";
    }

    static ParametersDescription description() {
      auto desc = EventImporter::description();
      desc.setDescription("HepMC3 LHEF file importer module");
      desc.add<std::string>("filename", "input.lhef").setDescription("Input filename");
      return desc;
    }

    bool operator>>(Event& evt) override {
      if (!reader_->readEvent())
        return false;
      evt.clear();
      const auto& hepeup = reader_->hepeup;
      CG_DEBUG("LHEFHepMCImporter:next").log([&hepeup](auto& log) { hepeup.print(log.stream()); });
      int id_ip1 = -1, id_ip2 = -1;
      pdgid_t pdg_ip1{0}, pdg_ip2{0};
      for (int i = 0; i < hepeup.NUP; ++i) {
        Particle part;
        part.setRole(Particle::Role::CentralSystem);
        part.setPdgId((long)hepeup.IDUP.at(i));
        const auto& hepeup_mom = hepeup.PUP.at(i);
        part.setMomentum(Momentum::fromPxPyPzE(hepeup_mom.at(0), hepeup_mom.at(1), hepeup_mom.at(2), hepeup_mom.at(3)),
                         false);
        part.setStatus(hepeup.ISTUP.at(i) < 0 ? Particle::Status::Propagator : Particle::Status::FinalState);
        if (hepeup.ISTUP.at(i) == -9) {
          part.setStatus(Particle::Status::PrimordialIncoming);
          if (part.momentum().pz() > 0) {
            part.setRole(Particle::Role::IncomingBeam1);
            id_ip1 = i;
            pdg_ip1 = hepeup.IDUP.at(i);
          } else {
            part.setRole(Particle::Role::IncomingBeam2);
            id_ip2 = i;
            pdg_ip2 = hepeup.IDUP.at(i);
          }
        }
        const auto& moth = hepeup.MOTHUP.at(i);
        if (moth.first > 0)
          part.addMother(evt[moth.first - 1]);
        if (moth.second > 0)
          part.addMother(evt[moth.second - 1]);
        if (utils::contains(part.mothers(), id_ip1)) {
          if (evt[Particle::Role::OutgoingBeam1].empty() && hepeup.IDUP.at(i) == (long)pdg_ip1)
            part.setRole(Particle::Role::OutgoingBeam1);
          else
            part.setRole(Particle::Role::Parton1);
        }
        if (utils::contains(part.mothers(), id_ip2)) {
          if (evt[Particle::Role::OutgoingBeam2].empty() && hepeup.IDUP.at(i) == (long)pdg_ip2)
            part.setRole(Particle::Role::OutgoingBeam2);
          else
            part.setRole(Particle::Role::Parton2);
        }
        evt.addParticle(part);
      }
      return true;
    }

  private:
    void initialise() override {
      const auto& heprup = reader_->heprup;
      CG_DEBUG("LHEFHepMCImporter").log([&heprup](auto& log) { heprup.print(log.stream()); });
      setCrossSection(Value{heprup.XSECUP.at(0), heprup.XERRUP.at(0)});
    }

    std::unique_ptr<LHEF::Reader> reader_;
  };
}  // namespace cepgen
REGISTER_EVENT_IMPORTER("lhef_hepmc", LHEFHepMCImporter);
