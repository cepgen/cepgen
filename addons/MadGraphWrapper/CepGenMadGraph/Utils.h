/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2024  Laurent Forthomme
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

#ifndef CepGenMadGraph_Utils_h
#define CepGenMadGraph_Utils_h

#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen::mg5amc {
  using ProcessParticles = std::pair<std::vector<std::string>, std::vector<std::string> >;
  /// Unpack the particles' content and role in the process from a string
  /// \param[in] proc string, human-built process definition
  ProcessParticles unpackProcessParticles(const std::string& proc);
  /// Unpack all particle properties from MadGraph
  /// \param[in] part_name mg5_aMC particle name
  /// \param[in] model mg5_aMC model to use
  ParticleProperties describeParticle(const std::string& part_name, const std::string& model = "");
  /// Run a mg5_aMC command and return its result
  /// \param[in] cmds list of commands to send to the mg5_aMC path
  /// \param[in] card_path filename to use for the steering card
  /// \param[in] keep_output keep the steering card after run?
  /// \return full mg5_aMC output
  std::vector<std::string> runCommand(const std::vector<std::string>& cmds,
                                      const std::string& card_path,
                                      bool keep_output = false);
  /// Normalise a process name to make it computer-readable
  /// \param[in] model mg5_aMC model to use
  std::string normalise(const std::string& proc, const std::string& model = "");
}  // namespace cepgen::mg5amc

#endif
