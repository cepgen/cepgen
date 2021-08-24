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

#ifndef CepGenAddOns_MadGraphWrapper_MadGraphInterface_h
#define CepGenAddOns_MadGraphWrapper_MadGraphInterface_h

#include <memory>
#include <string>

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/Filesystem.h"

// forward-declaration of base MadGraph standalone_cpp process
class CPPProcess;

namespace cepgen {
  class MadGraphInterface {
  public:
    MadGraphInterface(const ParametersList&);

    std::string run() const;

  private:
    static constexpr size_t cmd_buffer_size_ = 256;
    static const std::unordered_map<std::string, pdgid_t> mg5_parts_;

    static std::string runCommand(const std::string&);
    static std::string generateLibrary(const fs::path&, const fs::path&, const std::string&);
    static std::string generateProcess(const fs::path&);

    using ProcessParticles = std::pair<std::vector<pdgid_t>, std::vector<pdgid_t> >;
    static ProcessParticles unpackProcessParticles(const std::string&);

    void prepareCard() const;
    void linkCards() const;
    std::string prepareMadGraphProcess() const;

    const std::string proc_;
    const std::string model_;
    const fs::path card_path_;
    const fs::path standalone_cpp_path_;
    const fs::path tmp_dir_;
    const std::string log_filename_;
  };
}  // namespace cepgen

#endif
