/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

#ifndef CepGenMadGraph_Interface_h
#define CepGenMadGraph_Interface_h

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/Filesystem.h"

namespace cepgen::mg5amc {
  class Interface final : public SteeredObject<Interface> {
  public:
    explicit Interface(const ParametersList&);

    static ParametersDescription description();

    std::string run() const;

    /// Retrieve a CepGen-compatible parameters list from a MadGraph parameters card
    static ParametersDescription extractParamCardParameters(const std::string&);
    /// Generate a MadGraph parameters card from CepGen user-steered parameters
    static std::string generateParamCard(const ParametersDescription&);

  private:
    static constexpr size_t cmd_buffer_size_ = 256;
    static std::unordered_map<std::string, spdgid_t> mg5_parts_;

    static void generateLibrary(const fs::path&, const fs::path&, const fs::path&);

    void parseExtraParticles();
    void linkCards() const;
    std::string prepareMadGraphProcess() const;

    const std::string proc_;
    const std::string model_;
    const fs::path tmp_dir_;
    const fs::path card_path_;
    const fs::path log_filename_;
    const fs::path standalone_cpp_path_;
    const ParametersList extra_particles_, model_parameters_;

    std::string extra_part_definitions_;
  };
}  // namespace cepgen::mg5amc

#endif
